<!-- dx-header -->
# dnanexus_vcf_QC (DNAnexus Platform App)

## What does this app do?
This app runs a python script to perform simple vcf QC.

## What are typical use cases for this app?
This app should be executed stand-alone or as part of a DNAnexus workflow for a single sample.

## What data are required for this app to run?
The app requires a VCF file (.vcf) containing variants to be evaluated and a bed file defining the regions within which variants in the vcf should be evaluated.

## What does this app output?
The app outputs one file, where [outPrefix] is the vcf filename without extension:
1. [outPrefix].vcf.QC is a tab delimited file containing:
 - sample id
 - mean het ratio (mean AAF of het variants)
 - mean homo ratio (mean AAF of hom variants)
 - het:homo ratio (ratio of het to hom variants on autosomes)
 - X homo:het ratio (ratio of het to hom variants on chrX)
 - gender (inferred from X homo:het ratio)

## How does this app work?
The app runs a bash script which runs the python script which generated the output file. The output file is then uploaded to dnanexus.

## What are the limitations of this app
- Inferred gender is tuned to work with TSOE data and the TSOnePlus.bed. Use of other assays/bed files may result in incorrect inferred gender.
- Inferred gender may generate incorrect results for biological reasons e.g. for patient's whose parents are closely related, or patients with sex chromosome anomolies (Turner, Klinefelter etc.)

## This app was made by EMEE GLH
