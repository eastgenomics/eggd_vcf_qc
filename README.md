<!-- dx-header -->
# eggd_vcf_qc (DNAnexus Platform App)

## What does this app do?
This app runs a python script to calculate the mean het, mean hom and het:hom ratio of non reference variants in a given vcf.

## What are typical use cases for this app?
This app should be executed stand-alone or as part of a DNAnexus workflow for a single sample.

## What data are required for this app to run?
The app requires a VCF file (.vcf) containing variants to be evaluated and a bed file defining the regions within which variants in the vcf should be evaluated.

## What does this app output?
The app outputs one tab delimited file `[vcf_prefix].vcf.qc`, where `vcf_prefix` is the vcf filename without extension:

* sample id
* mean het ratio (mean AAF of het variants)
* mean hom ratio (mean AAF of hom variants)
* het:hom ratio (ratio of het to hom variants on autosomes)
* X het:hom ratio (ratio of het to hom variants on chrX)


## This app was made by EMEE GLH

