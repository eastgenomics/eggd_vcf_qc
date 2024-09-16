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
        'x_hom': [],
        'crap_homs': [],
        'crap_hets': []
    }

    for record in vcf_handle.fetch():
        sample_fields = record.samples[sample]

        if not all(x in sample_fields for x in ['AD', 'DP', 'GT']):
            # TODO - do we still want to do this? does this even happen?
            print(f"Missing field(s)")
            continue

        # using the sum of all allele depths instead of the format AD
        # field to be the informative read depths supporting each allele
        informative_total_depth = sum(sample_fields['AD'])
        non_ref_depth = sum(sample_fields['AD'][1:])

        non_ref_aaf = round(non_ref_depth / informative_total_depth, 4)

        if len(set(sample_fields['GT'])) == 1:
            # homozygous variant
            counts['hom'].append(non_ref_aaf)

            if sample_fields['GQ'] < 10:
                counts['crap_homs'].append(sample_fields['GQ'])

            if re.match(r'(chr)?x', record.chrom, re.IGNORECASE):
                counts['x_hom'].append(non_ref_aaf)
        else:
            # heterozygous variant
            counts['het'].append(non_ref_aaf)

            if sample_fields['GQ'] < 30:
                counts['crap_hets'].append(sample_fields['GQ'])

            if re.match(r'(chr)?x', record.chrom, re.IGNORECASE):
                counts['x_het'].append(non_ref_aaf)

        # handy print for the logs for sense checking
        print(
            f"GT: {sample_fields['GT']}\tADs: {sample_fields['AD']}\t\t"
            f"AD DP: {sum(sample_fields['AD'])}\tFMT_DP: {sample_fields['DP']}"
            f"\tAAF: {non_ref_aaf}"
        )

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
    if not counts['het'] and counts['hom']:
        # we don't have both het and hom variants => TODO figure out what to do
        return

    ratios = {}

    ratios['mean_het'] = sum(counts['het']) / len(counts['het'])
    ratios['mean_hom'] = sum(counts['hom']) / len(counts['hom'])

    ratios['het_hom_ratio'] = len(counts['het']) / len(counts['hom'])

    ratios['x_het_hom_ratio'] = None

    if counts['x_hom'] and counts['x_het']:
        ratios['x_het_hom_ratio'] = len(counts['x_het']) / len(counts['x_hom'])

    print(
        f"\nTotal het counts: {len(counts['het'])}"
        f"Total hom counts: {len(counts['hom'])}"
        f"Total x het counts: {len(counts['x_het'])}"
        f"Total x hom counts: {len(counts['x_hom'])}"
        f"Crap homs: {len(counts['crap_homs'])}"
        f"Crap hets: {len(counts['crap_hets'])}\n"
    )

    for field, value in ratios.items():
        print(f"{field}\t{value}")

    return ratios


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
        for k, v in ratios.items():
            fh.write(f'{k}\t{v}')


def upload_output_file(outfile) -> None:
    """
    Upload output file to project

    Parameters
    ----------
    outfile : str
        name of file to upload
    """
    url_file = dxpy.upload_local_file(
        outfile,
        folder=dxpy.bindings.dxjob.DXJob(
            os.environ.get('DX_JOB_ID')).describe()['folder']
    )

    return {"output_file": dxpy.dxlink(url_file)}


@dxpy.entry_point('main')
def main(vcf_file, bed_file, outfile=None):

    tmp_vcf = intersect_vcf_with_bed(vcf=vcf_file, bed=bed_file)
    het_hom_counts = get_het_hom_counts(tmp_vcf)
    ratios = calculate_ratios(het_hom_counts)

    if not outfile:
        outfile = f"{re.sub(r'.vcf(.gz)?$', '', vcf_file)}.vcf.qc"

    write_output_file(outfile=outfile, ratios=ratios)
    uploaded_file = upload_output_file(outfile=outfile)

    return uploaded_file


if os.path.exists('/home/dnanexus'):
    dxpy.run()
else:
    main(vcf_file=sys.argv[1], bed_file=sys.argv[2])
