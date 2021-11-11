#!/usr/bin/env python3

# Purpose: Reformat the default VCF output by Varscan2. This script will update the QUAl for each variant to be equal to the GQ (Genotype Quality).
# Authors: Max Marin (mgmarin@g.harvard.edu),

# Usage Example: Varscan.VCF.Updater.GQtoQUALpy.py --input M0011368_9.IllPE.bwa-mem.Varscan2.vcf > M0011368_9.IllPE.bwa-mem.Varscan2.QualUpdated.vcf



__doc__ = """Reformat the default VCF output by Varscan2. This script will update the QUAl for each variant to be equal to the GQ (Genotype Quality)."""

import sys
import argparse
import pysam


def main():

    #### Setup argparser
    parser = argparse.ArgumentParser(
        description="Update a Varscan2 VCF to have the GQ field as its QUAL tag ")

    parser.add_argument('--input', type=str, 
                        help="Path to input Varscan2 VCF")

    args = parser.parse_args()


    #### Set input and output PATHs ####

    i_VS2_VCF_PATH = args.input


    #### Update QUAL score to Genotype Quality (GQ) assigned by Varscan2 ####

    # Make pysam VCF object for input
    input_VCF = pysam.VariantFile(i_VS2_VCF_PATH, "r")


    # Print the header to STDOUT
    print(input_VCF.header, end='')


    # Set QUAL score to GQ field for each variant, then output variant

    for variant in input_VCF:    

        i_GQ = variant.samples[0]['GQ']

        variant.qual = i_GQ
        
        # Each updated variant entry is printed to standard-out
        print(variant, end='')



    #########################################


if __name__ == "__main__":
    sys.exit(main())
