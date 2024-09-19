from glob import glob
import os
import re
import subprocess
import sys

if os.path.exists('/home/dnanexus'):
    # running in DNAnexus
    subprocess.check_call([
        'pip', 'install', "--no-index", "--no-deps"
    ] + glob("packages/*"))

import dxpy
import pysam


def intersect_vcf_with_bed(vcf, bed) -> str:
    """
    Intersect the vcf with the given bed file, outputting the unique
    variants present in vcf

    Parameters
    ----------
    vcf : str
        vcf file to intersect with
    bed : str
        bed file to intersect against

    Returns
    -------
    str
        file name of intersected vcf

    Raises
    ------
    AssertionError
        Raised when non zero exit code returned from bedtools intersect
    """
    tmp_vcf = f"{vcf}.tmp"
    command = (
        f"bedtools intersect -header -wa -u -a {vcf} -b {bed} > {vcf}.tmp"
    )

    process = subprocess.run(
        command,
        shell=True,
        check=True,
        capture_output=True
    )

    assert process.returncode == 0, (
        f"Error in calling bedtools intersect: {process.stderr.decode()}"
    )

    return tmp_vcf


def get_het_hom_counts(vcf) -> dict:
    """
    Get the AD counts for het and hom variants from vcf file

    Parameters
    ----------
    vcf : str
        filename of vcf to get counts from

    Returns
    -------
    dict
        dict of per variant AAF split by het and hom

    Raises
    ------
    AssertionError
        Raised if more than one sample present in vcf
    """
    vcf_handle = pysam.VariantFile(vcf)
    samples = vcf_handle.header.samples

    # ensure we don't get passed a multisample vcf and therefore we can
    # just zero index sample fields from here on
    assert len(samples) == 1, "more than one sample present in vcf"

    sample = samples[0]

    counts = {
        'het': [],
        'hom': [],
        'x_het': [],
        'x_hom': []
    }

    for record in vcf_handle.fetch():
        sample_fields = record.samples[sample]

        if not all(x in sample_fields for x in ['AD', 'DP', 'GT']):
            # TODO - do we still want to do this? does this even happen?
            print(f"Missing field(s)")
            continue

        # TODO - should we skip ref sites incase a gvcf is provided?
        if sample_fields['GT'] == (0,0):
            continue

        # using the sum of all allele depths instead of the format AD
        # field to be the informative read depths supporting each allele
        informative_total_depth = sum(sample_fields['AD'])
        non_ref_depth = sum(sample_fields['AD'][1:])

        non_ref_aaf = non_ref_depth / informative_total_depth

        if len(set(sample_fields['GT'])) == 1:
            # homozygous variant
            counts['hom'].append(non_ref_aaf)

            if re.match(r'(chr)?x', record.chrom, re.IGNORECASE):
                counts['x_hom'].append(non_ref_aaf)
        else:
            # heterozygous variant
            counts['het'].append(non_ref_aaf)

            if re.match(r'(chr)?x', record.chrom, re.IGNORECASE):
                counts['x_het'].append(non_ref_aaf)

        # handy print for the logs for sense checking
        print(
            f"GT: {sample_fields['GT']}\tADs: {sample_fields['AD']}\t\t"
            f"AD DP: {sum(sample_fields['AD'])}\tFMT_DP: {sample_fields['DP']}"
            f"\tAAF: {non_ref_aaf}"
        )

    counts = {k: sorted(v) for k, v in counts.items()}

    return counts


def calculate_ratios(counts) -> dict:
    """
    Calculate the het hom ratios from the given counts

    Parameters
    ----------
    counts : dict
        dict of AAFs for het, hom, x_het and x_hom variants

    Returns
    -------
    dict
        dict containing calculated het hom ratios
    """
    ratios = {
        'mean_het': None,
        'mean_hom': None,
        'het_hom_ratio': None,
        'x_het_hom_ratio': None
    }

    if not counts['het'] or not counts['hom']:
        # we don't have both het and hom variants => TODO figure out what to do
        return ratios

    ratios['mean_het'] = f"{sum(counts['het']) / len(counts['het']):.4f}"
    ratios['mean_hom'] = f"{sum(counts['hom']) / len(counts['hom']):.4f}"

    ratios['het_hom_ratio'] = f"{len(counts['het']) / len(counts['hom']):.4f}"

    if counts['x_hom'] and counts['x_het']:
        ratios['x_het_hom_ratio'] = (
            f"{len(counts['x_het']) / len(counts['x_hom']):.4f}"
        )

    print(
        f"\nTotal het counts: {len(counts['het'])}\n"
        f"Total hom counts: {len(counts['hom'])}\n"
        f"Total x het counts: {len(counts['x_het'])}\n"
        f"Total x hom counts: {len(counts['x_hom'])}\n"
    )

    for field, value in ratios.items():
        print(f"{field}\t{value}")

    return ratios


def download_input_file(remote_file) -> str:
    """
    Download given input file with same name as file in project

    Parameters
    ----------
    remote_file : dict
        DNAnexus input file

    Returns
    -------
    str
        name of locally downloaded file
    """
    local_name = dxpy.describe(remote_file).get('name')
    dxpy.bindings.dxfile_functions.download_dxfile(
        dxid=remote_file,
        filename=local_name
    )

    return local_name


def write_output_file(outfile, ratios) -> None:
    """
    Write ratios to output file

    Parameters
    ----------
    outfile : str
        name of file to write to
    ratios : dict
        dict of field and calculated values to write
    """
    with open(outfile, 'w') as fh:
        print(ratios)
        header = '\t'.join(ratios.keys())
        values = '\t'.join([str(x) for x in ratios.values()])

        fh.write(f"{header}\n{values}\n")


def upload_output_file(outfile) -> None:
    """
    Upload output file to set folder in current project

    Parameters
    ----------
    outfile : str
        name of file to upload
    """
    print(f"Uploading {outfile}")

    url_file = dxpy.upload_local_file(
        filename=outfile,
        folder=dxpy.bindings.dxjob.DXJob(
            os.environ.get('DX_JOB_ID')).describe()['folder'],
        wait_on_close=True
    )

    return {"output_file": dxpy.dxlink(url_file)}


@dxpy.entry_point('main')
def main(vcf_file, bed_file):

    local_vcf_file = download_input_file(vcf_file)
    local_bed_file = download_input_file(bed_file)

    tmp_vcf = intersect_vcf_with_bed(vcf=local_vcf_file, bed=local_bed_file)
    het_hom_counts = get_het_hom_counts(tmp_vcf)
    ratios = calculate_ratios(het_hom_counts)

    outfile = f"{re.sub(r'.vcf(.gz)?$', '', local_vcf_file)}.vcf.qc"

    if os.path.exists('/home/dnanexus'):
        write_output_file(outfile=outfile, ratios=ratios)
        uploaded_file = upload_output_file(outfile=outfile)

        return uploaded_file


if os.path.exists('/home/dnanexus'):
    dxpy.run()
elif __name__ == "__main__":
    main(vcf_file=sys.argv[1], bed_file=sys.argv[2])
