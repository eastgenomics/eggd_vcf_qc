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


"""
TODO
- intersect vcf with bed and get unique variants
- calculate het hom ratio
- calculate hom het ratio

- loop over records:
  - loop over samples:
    - skip samples with no AD or DP or GT
    - get ref and alt (consider multi allelics here)
    - consider indels not having ref or alt depths (app sets count and dp to 1)
    - AAF calculated as alt_count*1.0/int(record.samples[sample]['DP'])
    - if len GT == 1:
        - add AAD as HOM
    - else:
        - add AAF as HET

"""


def intersect_vcf_with_bed(vcf, bed) -> str:
    """
    Intersect the vcf with the given bed file

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

    process = subprocess.run(command, shell=True, capture_output=True)

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
    str
        string of sample name parsed from vcf header
    dict
        dict of per variant AAF split by het and hom

    Raises
    ------
    AssertionError
        Raised when more than one sample present in vcf
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

        # using the sum of all allele depths instead of the format AD
        # field to be the informative read depths supporting each allele
        informative_total_depth = sum(sample_fields['AD'])
        non_ref_depth = sum(sample_fields['AD'][1:])

        non_ref_aaf = round(non_ref_depth / informative_total_depth, 4)

        print(
            f"GT: {sample_fields['GT']}\tADs: {sample_fields['AD']}\t\t"
            f"AD DP: {sum(sample_fields['AD'])}\tFMT_DP: {sample_fields['DP']}"
            f"\tAAF: {non_ref_aaf}"
        )

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

    return sample, counts


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

    print(f"\nTotal het counts: {len(counts['het'])}")
    print(f"Total hom counts: {len(counts['hom'])}")
    print(f"Total x het counts: {len(counts['x_het'])}")
    print(f"Total x hom counts: {len(counts['x_hom'])}\n")

    return ratios



@dxpy.entry_point('main')
def main(vcf_file, bed_file, outfile=None):

    tmp_vcf = intersect_vcf_with_bed(vcf=vcf_file, bed=bed_file)
    sample, het_hom_counts = get_het_hom_counts(tmp_vcf)
    ratios = calculate_ratios(het_hom_counts)

    for x, y in ratios.items():
        print(x, "\t", y)


if os.path.exists('/home/dnanexus'):
    dxpy.run()

if __name__ == "__main__":
    main(vcf_file=sys.argv[1], bed_file=sys.argv[2])
