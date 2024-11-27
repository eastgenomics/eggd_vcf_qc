<!-- dx-header -->
# eggd_vcf_qc (DNAnexus Platform App)

## What does this app do?
This app runs a python script to calculate the mean het, mean hom and het:hom ratio of non reference variants in a given VCF. Alternate heterozygous variants (i.e. those with a genotype of `1/2`) are not included in the mean het or het:hom ratios.

## What are typical use cases for this app?
This app should be executed stand-alone or as part of a DNAnexus workflow for a single sample.

## What inputs are required for this app to run?
##### Required
- `-ivcf_file` (`file`): VCF file containing variants to be evaluated.
- `-ibed_file` (`file`): BED file defining the regions within which variants in the VCF should be evaluated.

### Example commands
Running on DNAnexus:
```
dx run eggd_vcf_qc \
    -ivcf_file=file-xxx \
    -ibed_file=file-xxx
```

Running locally:
```
python3 src/vcf_qc.py <vcf-file> <bed-file>
```

## What does this app output?
- `output_file` (`file`): Tab-delimited file  `[vcf_prefix].vcf.QC`, where `vcf_prefix` is the VCF filename without extension. This includes columns:

Column | Description |
--- | --- |
Sample | Sample ID |
mean_het | Mean AAF of het variants |
mean_hom | Mean AAF of hom variants |
het_hom_ratio | Ratio of het to hom variants on autosomes
x_het_hom_ratio | Ratio of het to hom variants on chrX |

## This app was made by EMEE GLH

