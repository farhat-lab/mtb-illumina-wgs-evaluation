#!/usr/bin/env python3

# Purpose: Script for parsing and converting NucDiff SV GFF format to a custom BED format, along with processing and filtering.
# Authors: Max Marin (mgmarin@g.harvard.edu),



__doc__ = """Script for parsing and converting NucDiff SV GFF format to a custom BED format, along with processing and filtering"""

import sys
import os
import argparse
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Parse and merge relevant sample RNAseq expression,"
                    "quality control and metadata for a given experiment "
                    "(FEX) and gene cluster")
    parser.add_argument('--nucdiff_sv_gff', type=str, 
                        help="Path to input NucDiff_Filtered_GFF_PATH")


    parser.add_argument('--sample_id', type=str,
                        help="Sample ID to add to output table")

    args = parser.parse_args()



    #### Set input parameters and PATHs ####

    # A) Define Experiment


    i_NucDiff_Query_Ref_Struct_GFF_PATH = args.nucdiff_sv_gff

    sampleID = args.sample_id

    #########################################



    # Define output PATHs 

    NucDiff_Filtered_SVs_Output_Prefix = i_NucDiff_Query_Ref_Struct_GFF_PATH.rstrip(".gff") 

    i_NucDiff_Filtered_SVs_BED_PATH = NucDiff_Filtered_SVs_Output_Prefix + ".Parsed.bed"  

    i_NucDiff_Filtered_SVs_BED_With50bpLengthFilter_PATH = NucDiff_Filtered_SVs_Output_Prefix + ".Parsed.With50bpLengthFilter.bed"

    i_NucDiff_Filtered_SVs_BED_With50bpLengthFilter_NoInversions_PATH = NucDiff_Filtered_SVs_Output_Prefix + ".Parsed.With50bpLengthFilter.NoInversions.bed"




    i_NucDiff_SV_DF = parse_NucDiff_SV_GFF(i_NucDiff_Query_Ref_Struct_GFF_PATH)
    
    i_NucDiff_SV_DF.columns = ["Chrom", "start_0based", "end_0based", "SV_Type", "SV_Length"]

    i_NucDiff_SV_DF["SampleID"] = sampleID
    
    i_NucDiff_SV_DF = i_NucDiff_SV_DF.sort_values(["start_0based", "end_0based"])
    
    i_NucDiff_SV_Wi50bp_LenFilter_DF = i_NucDiff_SV_DF.query("SV_Length >= 50")
    
    i_NucDiff_SV_Wi50bp_LenFilter_NoInversions_DF = i_NucDiff_SV_Wi50bp_LenFilter_DF[ i_NucDiff_SV_Wi50bp_LenFilter_DF["SV_Type"] != "inversion" ] 
    



    i_NucDiff_SV_DF.to_csv(i_NucDiff_Filtered_SVs_BED_PATH, sep = "\t", index = False, header = False)
    
    i_NucDiff_SV_Wi50bp_LenFilter_DF.to_csv(i_NucDiff_Filtered_SVs_BED_With50bpLengthFilter_PATH, sep = "\t", index = False, header = False)

    i_NucDiff_SV_Wi50bp_LenFilter_NoInversions_DF.to_csv(i_NucDiff_Filtered_SVs_BED_With50bpLengthFilter_NoInversions_PATH, sep = "\t", index = False, header = False)
    







def GFF_AttributesCol_To_Dict(input_GFF_AttributesCol_Str):
    
    ###For each line of your VCF, create a dictionnary with this array key : info, value : value of this info
    dict_INFO = {}
    for i in input_GFF_AttributesCol_Str.split(";"):
        ###Just looking for line with "=" character (as key = value)
        if "=" in i:
            ###Left from equal sign is key (Gene.refGene, ExonicFunc.refGene...)
            key = i.split("=")[0]
            ###Right from equal sign is value (RBL1,synonymous_SNV...)
            value = i.split("=")[1]
            ###Put them in a dictionnary
            dict_INFO[key]=value
  
    return dict_INFO
        


def parse_NucDiff_SV_GFF(input_NucDiff_Query_Ref_Struct_GFF_PATH):

    SV_CNT = 1

    All_Valid_SV_Types = ['insertion', 'duplication', 'tandem_duplication', 'deletion',  'collapsed_repeat',  'collapsed_tandem_repeat', 'inversion', 'substitution']
    INS_SV_Types = ['insertion', 'duplication',       'tandem_duplication']
    DEL_SV_Types = ['deletion',  'collapsed_repeat',  'collapsed_tandem_repeat']
    SUB_SV_Types = ['substitution']
    BLK_SV_Types = ['inversion']

    listOf_InfoTuples = []


    with open(input_NucDiff_Query_Ref_Struct_GFF_PATH) as input_SV_GFF:

        for line in input_SV_GFF:
            if not line.startswith("#"):
                NucDiff_GFF_Row_Line = line.rstrip("\n").split("\t")

                i_seqID_Chrom = NucDiff_GFF_Row_Line[0]
                i_source = NucDiff_GFF_Row_Line[1]
                i_type = NucDiff_GFF_Row_Line[2]            
                i_start = int( NucDiff_GFF_Row_Line[3] )    
                i_end = int( NucDiff_GFF_Row_Line[4] )
                i_score = NucDiff_GFF_Row_Line[5]
                i_strand = NucDiff_GFF_Row_Line[6]
                i_phase = NucDiff_GFF_Row_Line[7]
                i_attributes = NucDiff_GFF_Row_Line[8]

                i_Attributes_Dict = GFF_AttributesCol_To_Dict(i_attributes)

                NucDiff_SV_Type = i_Attributes_Dict["Name"]

                if not NucDiff_SV_Type in All_Valid_SV_Types: continue


                if "ins_len"  in i_Attributes_Dict.keys():
                    GEN_SV_Len = i_Attributes_Dict["ins_len"]
                    
                    
                    

                elif "del_len"  in i_Attributes_Dict.keys():
                    GEN_SV_Len = i_Attributes_Dict["del_len"]

                elif "subst_len"  in i_Attributes_Dict.keys():
                    GEN_SV_Len = i_Attributes_Dict["subst_len"]

                elif "blk_len"  in i_Attributes_Dict.keys():
                    GEN_SV_Len = i_Attributes_Dict["blk_len"]

                else: 
                    GEN_SV_Len = 0
                    print(NucDiff_SV_Type, "has no length attribute")


                GEN_SV_Len = int(GEN_SV_Len)
                
                
                if NucDiff_SV_Type in ['duplication', 'tandem_duplication']: # Shift start of all duplications to the left by the length of the duplication (For NucDiff Duplications)
                    i_start = i_start - GEN_SV_Len
                

                
                i_start_0based = i_start - 1
                i_end_0based = i_end
                
                SV_Info_Tuple = (i_seqID_Chrom, i_start_0based, i_end_0based, NucDiff_SV_Type, GEN_SV_Len, )


                #if (not "ins_len"  in i_Attributes_Dict.keys()) and (not "del_len"  in i_Attributes_Dict.keys()) :
                #if GEN_SV_Len >= 5000:
                #    print(SV_Info_Tuple)    
                #    print(NucDiff_GFF_Row_Line)
                #    print("")

                SV_CNT += 1

                #if SV_CNT >= 15: break

                listOf_InfoTuples.append(SV_Info_Tuple)


    NucDiff_SV_Info_DF = pd.DataFrame(listOf_InfoTuples)

    return NucDiff_SV_Info_DF




if __name__ == "__main__":
    sys.exit(main())

