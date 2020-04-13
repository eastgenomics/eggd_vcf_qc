#!/usr/bin/python2
# 
# 
# 
# 
# Kim Brugger (21 Mar 2018), contact: kim@brugger.dk
# Matt Garner (20191107)

import sys
import subprocess
import os
import argparse
import pysam

def analyse_vcf( vcf_file, bed_file ):
    """ Does the QC of a vcf file

    Args:
      vcf_file(str): name of vcf file

    Returns:
      None
    """

    # Make temp on-target vcf
    temp_vcf = vcf_file + ".tmp"
    command = "bedtools intersect -header -a {vcf_file} -b {bed_file} -wa > {temp_vcf}".format(vcf_file=vcf_file, bed_file=bed_file, temp_vcf=temp_vcf)

    subprocess.call(command, shell=True)
    
    vcf = pysam.VariantFile( temp_vcf )
    samples = vcf.header.samples

    X_splits = {}
    X_splits['HOMO'] = []
    X_splits['HET' ] = []

    Y_splits = {}
    Y_splits['HOMO'] = []
    Y_splits['HET' ] = []
  
    homo_splits = []
    het_splits  = []


    splits = {}
    for sample in samples:
        splits[ sample ] = {}
        splits[ sample ][ 'HET'] = []
        splits[ sample ][ 'HOM'] = []
        splits[ sample ][ 'X' ]  = {}
        splits[ sample ][ 'X' ][ 'HET'] = []
        splits[ sample ][ 'X' ][ 'HOM'] = []

    for record in vcf.fetch():
        for sample in samples:

            if ( 'AD' not in record.samples[sample] or 
                 'DP' not in record.samples[sample] ):
                continue

            gts = record.samples[sample]['GT']
            counts = record.samples[sample]['AD']
            if not gts:
                continue
            if not counts:
                continue

            try:
                ref_count, alt_count = counts[ gts[0]], counts[ gts[1]]
            except:
                continue
            
            # Indels sometimes have no reference and alt depths (for some reason!)
            if ( record.samples[sample]['DP'] == 0 and ref_count == 0 and alt_count == 0):
                alt_count = 1
                record.samples[sample]['DP'] = 1

            AAF = alt_count*1.0/int(record.samples[sample]['DP'])

            if len(set(record.samples[sample]['GT'])) == 1:
                splits[ sample ][ 'HOM' ].append( AAF )

                if record.chrom == 'X':
                    splits[ sample ][ 'X' ][ 'HOM'].append( AAF )
            else:
                splits[ sample ][ 'HET' ].append( AAF )

                if record.chrom == 'X':
                    splits[ sample ][ 'X' ][ 'HET'].append( AAF )


    for sample in samples:

        if (len( splits[ sample ][ 'HET']) > 0) and (len( splits[ sample ][ 'HOM']) > 0) :
            mean_het = sum( splits[ sample ][ 'HET'])*1.0/len( splits[ sample ][ 'HET'] )
            mean_hom = sum( splits[ sample ][ 'HOM'])*1.0/len( splits[ sample ][ 'HOM'] )
            het_hom_ratio = len( splits[ sample ][ 'HET'] )*1.0/len( splits[ sample ][ 'HOM'] )

        # If no hets/homs to calculate ratios/scores, set values to -1000
        else:
            mean_het = -1000
            mean_hom = -1000
            het_hom_ratio = -1000

        if ( len(splits[ sample ][ 'X' ][ 'HOM']) == 0):
            X_het_hom_ratio = -1000
        else:
            X_het_hom_ratio = len(splits[ sample ][ 'X' ][ 'HET'])*1.0/len(splits[ sample ][ 'X' ][ 'HOM'])

        if (X_het_hom_ratio > 1 ):
            gender = "male"
        else:
            gender = "female"

        mean_het = "{:.4f}".format( mean_het )
        mean_hom = "{:.4f}".format( mean_hom )
        het_hom_ratio = "{:.4f}".format( het_hom_ratio )
        X_het_hom_ratio = "{:.4f}".format( X_het_hom_ratio )

        print( "\t".join( [ sample, mean_het, mean_hom, het_hom_ratio, X_het_hom_ratio, gender]))

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='vcf_QC.py: simple QC of a vcf file ')

    parser.add_argument('-o', '--outfile', help="if set writes to this file")
    parser.add_argument('vcf_files', metavar='vcf-files', nargs="+",  help="vcf file(s) to analyse")
    parser.add_argument('bed_file', metavar='bed-file', help="bed file containing regions to analyse")
    
    args = parser.parse_args()

    if args.outfile:
        sys.stdout = open("{}".format( args.outfile), 'w')

    print "\t".join([ "Sample", "mean het ratio", "mean homo ratio", "het:homo ratio", "X homo:het ratio", "gender"])

    for vcf_file in args.vcf_files:
        analyse_vcf( vcf_file, args.bed_file )