#!/usr/bin/python


### EBR Calculation Script - V7
### Maximillian Marin (mgmarin@g.harvard.edu)




## Import statements
import sys
import json
import pandas as pd
import numpy as np




# Run the script as follows: python _____.py <> <> <>

# arg1 = <GroundTruth.minimap2.mpileup.call.vcf>
# arg2 = <GroundTruth.NucDiff.Deletions.gff>
# arg3 = <Query.Illumina.Pilon.VCF>
# arg4 = <Output.EBR.IndividualIsolate.tsv>
# arg5 = <OutcomeClassifcations.json>




### Define functions for EBR calculation

def parse_PB_mpileup_ToDict(i_MM2_AtoRef_PP_BAM_mpileup_out_PATH):
    i_PB_mpileup_Dict = {}

    prev_POS = 0

    with open(i_MM2_AtoRef_PP_BAM_mpileup_out_PATH, "r") as sample_MM2_mpileup_vcf:

        for line in (sample_MM2_mpileup_vcf): #tqdm

            if not line.startswith("#"):

                VCF_Row_Line = line.split("\t")

                CHROM = VCF_Row_Line[0]

                POS = int(VCF_Row_Line[1])

                REF = VCF_Row_Line[3]

                ALT = VCF_Row_Line[4]

                QUAL = VCF_Row_Line[5]

                FILTER = VCF_Row_Line[6]

                INFO = VCF_Row_Line[7]
                #INFO_Split = INFO.split(";")

                ###For each line of your VCF, create a dictionnary with this array key : info, value : value of this info
                dict_INFO = {}
                for i in INFO.split(";"):
                    ###Just looking for line with "=" character (as key = value)
                    if "=" in i:
                        ###Left from equal sign is key (Gene.refGene, ExonicFunc.refGene...)
                        key = i.split("=")[0]
                        ###Right from equal sign is value (RBL1,synonymous_SNV...)
                        value = i.split("=")[1]
                        ###Put them in a dictionnary
                        dict_INFO[key]=value

                row_DP = int( dict_INFO["DP"] )

                if ALT == "<*>":
                    ALT = "."
                else:
                    ALT = ALT.rstrip(",<*>")

                #i_PB_mpileup_Dict[POS] = (POS, REF, ALT, row_DP)

                if POS != prev_POS:
                    i_PB_mpileup_Dict[POS] = (POS, REF, ALT, row_DP)

                prev_POS = POS
                
    return i_PB_mpileup_Dict


def parse_Ill_PilonVCF_ToDict(i_Ill_Pilon_VCF_PATH):

    i_Ill_Pilon_VCF_Dict = {}

    prev_POS = 0

    with open(i_Ill_Pilon_VCF_PATH, "r") as sample_Ill_Pilon_vcf:

        for line in (sample_Ill_Pilon_vcf): #tqdm

            if not line.startswith("#"):

                VCF_Row_Line = line.split("\t")

                CHROM = VCF_Row_Line[0]

                POS = int(VCF_Row_Line[1])

                REF = VCF_Row_Line[3]

                ALT = VCF_Row_Line[4]

                QUAL = VCF_Row_Line[5]

                FILTER = VCF_Row_Line[6]

                INFO = VCF_Row_Line[7]
                #INFO_Split = INFO.split(";")
                

                ###For each line of your VCF, create a dictionnary with this array key : info, value : value of this info
                dict_INFO = {}
                for i in INFO.split(";"):
                    ###Just looking for line with "=" character (as key = value)
                    if "=" in i:
                        ###Left from equal sign is key (Gene.refGene, ExonicFunc.refGene...)
                        key = i.split("=")[0]
                        ###Right from equal sign is value (RBL1,synonymous_SNV...)
                        value = i.split("=")[1]
                        ###Put them in a dictionnary
                        dict_INFO[key]=value

                row_DP = int( dict_INFO.get("DP", 0) )
                row_MQ = int( dict_INFO.get("MQ", 0) )
                row_TD = int( dict_INFO.get("TD", 0) )
                row_AF = float( dict_INFO.get("AF", 0) )                
                
                
                if POS != prev_POS:
                    i_Ill_Pilon_VCF_Dict[POS] = (POS, REF, ALT, FILTER, row_DP, row_MQ, row_TD, row_AF)
                    
                elif POS == prev_POS:
                    continue
                    #print("position occured twice!")
                    #print(i_Ill_Pilon_VCF_Dict[prev_POS] )
                    #print((POS, REF, ALT, FILTER, row_DP, row_MQ, row_TD, row_AF))

                prev_POS = POS

    return i_Ill_Pilon_VCF_Dict




def parse_NucDiff_SVs_DeletedPositions_FromGFF(NucDiff_SVs_Filtered_DeletedPositionsOnly_GFF_PATH):
    
    
    NucDiff_SVs_Dels_GFF_DF = pd.read_csv(NucDiff_SVs_Filtered_DeletedPositionsOnly_GFF_PATH, sep="\t", header=None)

    # NOTE: GFFs use 1-based coordinates
    NucDiff_SVs_Dels_GFF_DF.columns = ["Chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]

    # NucDiff_SVs_Dels_GFF_DF["Del_Length"] = NucDiff_SVs_Dels_GFF_DF["end"] - NucDiff_SVs_Dels_GFF_DF["start"] + 1

    i_PB_NucDiff_DeletionsDict = {}

    for i, row in NucDiff_SVs_Dels_GFF_DF.iterrows():
        #print(i, row)

        Del_Start_1based = row["start"]
        Del_End_1based = row["end"]

        #if Del_Start_1based == 698723: continue

        for ref_POS in np.arange(Del_Start_1based, Del_End_1based + 1):

            NucDiff_SV_Is_Deleted = True

            i_PB_NucDiff_DeletionsDict[ref_POS] = NucDiff_SV_Is_Deleted


    return i_PB_NucDiff_DeletionsDict
    




def get_AgreementTuple_Between_PB_and_Ill_V7(input_PB_mpileup_Dict, input_Ill_Pilon_VCF_Dict, input_PB_NucDiff_DeletionsDict):
    

    listOfTuples = []

    
    count_EBR_Outcomes_Dict = {}
    
    count_EBR_Outcomes_Dict["Outcome_1"] = 0
    count_EBR_Outcomes_Dict["Outcome_2"] = 0    
    count_EBR_Outcomes_Dict["Outcome_3"] = 0
    count_EBR_Outcomes_Dict["Outcome_4"] = 0   
    count_EBR_Outcomes_Dict["Outcome_5"] = 0
    count_EBR_Outcomes_Dict["Outcome_6"] = 0   
    count_EBR_Outcomes_Dict["Outcome_7"] = 0
    count_EBR_Outcomes_Dict["Outcome_8"] = 0    
    count_EBR_Outcomes_Dict["Outcome_9"] = 0    
    count_EBR_Outcomes_Dict["Outcome_10"] = 0    
    count_EBR_Outcomes_Dict["Outcome_11_Unknown"] = 0    

    

    count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_ValidCoverage_LowMQ"] = 0 
    count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_ValidCoverage_FlaggedPairs"] = 0 
    count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_TotalCoverage"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Amb"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Del"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Amb;LowCov"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;Amb"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;LowCov"] = 0
    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;Amb;LowCov"] = 0    
    count_EBR_Outcomes_Dict["Outcome_4_SO_Other"] = 0  
    
    
    
    

    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_ValidCoverage_LowMQ"] = 0 
    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_ValidCoverage_FlaggedPairs"] = 0 
    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_TotalCoverage"] = 0    
    count_EBR_Outcomes_Dict["Outcome_6_SO_Amb"] = 0    
    count_EBR_Outcomes_Dict["Outcome_6_SO_Del"] = 0    
    count_EBR_Outcomes_Dict["Outcome_6_SO_Amb;LowCov"] = 0    
    count_EBR_Outcomes_Dict["Outcome_6_SO_Del;Amb"] = 0    
    count_EBR_Outcomes_Dict["Outcome_6_SO_Del;LowCov"] = 0
    count_EBR_Outcomes_Dict["Outcome_6_SO_Del;Amb;LowCov"] = 0   
    count_EBR_Outcomes_Dict["Outcome_6_SO_Other"] = 0    


    
    
    
    count_EBR_Outcomes_Dict["Outcome_9_DP2_NtAgree"] = 0 
    count_EBR_Outcomes_Dict["Outcome_9_DP2_NtDisagree"] = 0 
        
    count_EBR_Outcomes_Dict["NAN_Counts"] = 0
    

    
    # Evaluate "Agreement"/EBR for every position of the H37rv reference
    
    
    for ref_POS in ( np.arange(1, 4411532 + 1, 1) ): #tqdm
        #print(ref_POS)

        ### STEP 1: Parse relevant info from PBMM2 and Pilon VCFs ###
        
        # Get PacBio-Minimap2-Mpileup-Call info from Dict
        PB_Pos, PB_REF, PB_ALT, PB_DP = input_PB_mpileup_Dict.get(ref_POS, (ref_POS, 0, 0, 0))

        if (PB_REF == 0) and (PB_ALT == 0):
            Within_PB_Alignment = False
        else:
            Within_PB_Alignment = True

        # Get Illumina-Pilon info from Dict
        Ill_Pos, Ill_REF, Ill_ALT, Ill_Filter, Ill_DP, Ill_MQ, Ill_TD, Ill_AF = input_Ill_Pilon_VCF_Dict[ref_POS]


        
        ### STEP 2: Parse DELETION status information from the NucDiff_DeletionsDict
        
        
        PB_NucDiff_SV_Is_Deleted = input_PB_NucDiff_DeletionsDict.get(ref_POS, False)
        
        
        # Determine whether analysis of the PacBio assembly supports that the position is deleted
        
        PB_MM2_Supports_Del = ((PB_DP == 0) and (Within_PB_Alignment == True))
        
        
        PB_Supports_Del = (PB_MM2_Supports_Del or PB_NucDiff_SV_Is_Deleted )
        
        
        
        
        ### STEP 3: Classify agreement outcome based on comparison of PB-MM2 and Pilon variant calling  ###
        
        ##### Potential Outcomes #####
        
        # 1) PB-MM2 and Pilon are both confident, and AGREE (EBR = 1)
            # Requirements: PB-MM2 DP = 1, Pilon-Tag = PASS, Genotypes Agree

        if ( (PB_DP == 1) and (Ill_Filter == "PASS") and (PB_REF == Ill_REF) and (PB_ALT == Ill_ALT) ) : # Both PB and Ill agree that basepair is present and agree for basecall    
            # Set EBR = 1
            listOfTuples.append( (ref_POS, 1, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["Outcome_1"] += 1             
    
    
    
    
    
        # 2) PB-MM2 and Pilon don't support presence of nucleotide,  (EBR = 1)
            # Requirements: PB-MM2 DP = 0, Pilon-Tag != PASS

            # Sub-outcomes:
                # a) Low-Cov - Insufficient Physical Coverage

                # b) Low-Cov - Insufficient reads passing heuristics

                # c) Amb

                # d) Del 
                
                # e) Multiple-Fail-Tags
       
        # Outcomes 2,3,4
        elif (PB_DP == 0)  and  (Ill_Filter != "PASS"): # Both PB and Ill agree that basepair is not present OR non-PASS
            
            
            listOf_PilonTags_WithDel = ["Del", "Del;LowCov", "Del;Amb", "Del;Amb;LowCov"]
            listOf_PilonTags_WithLowCov = ["LowCov", "Del;LowCov", "Amb;LowCov", "Del;Amb;LowCov"]
            
            ## Outcome 2) Both PacBio and Pilon agree that position is deleted 
            if (Ill_Filter in listOf_PilonTags_WithDel) and (PB_Supports_Del):

                listOfTuples.append( (ref_POS, 1, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
                count_EBR_Outcomes_Dict["Outcome_2"] += 1                  

            
            ## Outcome 3) PacBio supports deletion, but Pilon-Tag != (PASS or Del) 
            elif (Ill_Filter not in listOf_PilonTags_WithDel) and (PB_Supports_Del):
                
                listOfTuples.append( (ref_POS, 0, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
                count_EBR_Outcomes_Dict["Outcome_3"] += 1    
            
                #if (Ill_Filter in listOf_PilonTags_WithDel) and (PB_Supports_Del):
            
            
            ## Outcome 4) Both PacBio and Pilon do not support a given position, AND PacBio DOES NOT support a deletion at position
            ###

            # Sub-outcomes:
                # a) Low-Cov - Insufficient Physical Coverage

                # b) Low-Cov - Insufficient reads passing heuristics

                # c) Amb

                # d) Del 
                
                # e) Multiple-Fail-Tags


            elif (PB_Supports_Del == False):
                
                # Set EBR = np.nan (USED TO BE EBR = 0)
                listOfTuples.append( (ref_POS, np.nan, Ill_Filter))
                
                count_EBR_Outcomes_Dict["NAN_Counts"] += 1
                
                count_EBR_Outcomes_Dict["Outcome_4"] += 1  

                # Sub-outcomes

                if (Ill_Filter == "LowCov"):

                    if (Ill_TD >= 5) and (Ill_DP < 5) and (Ill_MQ == 0):
                        count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_ValidCoverage_LowMQ"] += 1

                    elif (Ill_TD >= 5) and (Ill_DP < 5) and (Ill_MQ > 0):
                        count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_ValidCoverage_FlaggedPairs"] += 1

                    elif (Ill_TD < 5) :
                        count_EBR_Outcomes_Dict["Outcome_4_SO_LowCov_Insufficient_TotalCoverage"] += 1 

                    else: print("Unknown LowCov outcome!!!!", ref_POS)          


                elif (Ill_Filter == "Amb"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Amb"] += 1 

                elif (Ill_Filter == "Del"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Del"] += 1 

                elif (Ill_Filter == "Amb;LowCov"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Amb;LowCov"] += 1 

                elif (Ill_Filter == "Del;Amb"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;Amb"] += 1 

                elif (Ill_Filter == "Del;LowCov"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;LowCov"] += 1 

                elif (Ill_Filter == "Del;Amb;LowCov"):
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Del;Amb;LowCov"] += 1 
                    
                elif (";" in Ill_Filter):
                    print("Outcome_4 - Not in dict, multiple tags", Ill_Filter)
                    print(input_Ill_Pilon_VCF_Dict[ref_POS])
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Other"] += 1 

                else:
                    print("Unkown sub-outcome for Outcome #4, agree that basepair is NOT confidently present")
                    count_EBR_Outcomes_Dict["Outcome_4_SO_Other"] += 1


        ## Outcome 5) PB-MM2 and Pilon are both confident, and DISAGREE (EBR = 0)
            # Requirements: PB-MM2 DP = 1, Pilon-Tag = PASS, Genotypes Disagree        

        elif ( (PB_DP == 1) and (Ill_Filter == "PASS") and (PB_ALT != Ill_ALT) ): # Both PB and Ill agree that basepair is present but DISAGREE for genotype    
            print("Genotype disagrement between PacBio and Illumina @", ref_POS)
            print(input_Ill_Pilon_VCF_Dict[ref_POS])
            print(input_PB_mpileup_Dict.get(ref_POS, (ref_POS, 0, 0, 0)))
            # Set EBR = 0
            listOfTuples.append( (ref_POS, 0, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["Outcome_5"] += 1             


                
        ## Outcome 6) PB-MM2 is confident, Pilon != PASS  (EBR = 0)
            # Requirements: PB-MM2 DP = 1, Pilon-Tag != PASS
            
            # Sub-outcomes:
                # a) Low-Cov - Insufficient Physical Coverage

                # b) Low-Cov - Insufficient reads passing heuristics

                # c) Amb

                # d) Del 
                
                # e) Multiple-Fail-Tags
        
        elif ( (PB_DP == 1) and (Ill_Filter != "PASS") ): #  
                        
            # Set EBR = 0
            listOfTuples.append( (ref_POS, 0, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["Outcome_6"] += 1             
  

            # Sub-outcomes
            
            if (Ill_Filter == "LowCov"):
                
                if (Ill_TD >= 5) and (Ill_DP < 5) and (Ill_MQ == 0):
                    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_ValidCoverage_LowMQ"] += 1

                elif (Ill_TD >= 5) and (Ill_DP < 5) and (Ill_MQ > 0):
                    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_ValidCoverage_FlaggedPairs"] += 1

                elif (Ill_TD < 5) :
                    count_EBR_Outcomes_Dict["Outcome_6_SO_LowCov_Insufficient_TotalCoverage"] += 1 

                else: print("Outcome_6 - Unknown LowCov outcome!!!!", ref_POS)          
                
                
                
            elif (Ill_Filter == "Amb"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Amb"] += 1 
                
            elif (Ill_Filter == "Del"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Del"] += 1 
                                
            elif (Ill_Filter == "Amb;LowCov"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Amb;LowCov"] += 1 

            elif (Ill_Filter == "Del;Amb"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Del;Amb"] += 1 
                
            elif (Ill_Filter == "Del;LowCov"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Del;LowCov"] += 1 
                
            elif (Ill_Filter == "Del;Amb;LowCov"):
                count_EBR_Outcomes_Dict["Outcome_6_SO_Del;Amb;LowCov"] += 1 
                
            elif (";" in Ill_Filter):
                print("Outcome_6 - Not in dict, multiple tags", Ill_Filter)
                count_EBR_Outcomes_Dict["Outcome_6_SO_Other"] += 1 


            else:
                print("Unkown sub-outcome for Outcome #6, agree that basepair is NOT confidently present")
                count_EBR_Outcomes_Dict["Outcome_6_SO_Other"] += 1
                

                

        ## Outcome 7) PB-MM2 DP = 0, PB supports deletion, but Pilon = PASS (EBR = 0)
            # Requirements: PB-MM2 DP = 0, Pilon-Tag == PASS, Within_PB_Alignment == True
            
        elif ( (PB_DP == 0) and (Ill_Filter == "PASS") and (PB_Supports_Del == True) ) :  

            # Set EBR = 1
            listOfTuples.append( (ref_POS, 0, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["Outcome_7"] += 1 



        ## Outcome 8) PB-MM2 DP = 0, PB DOES NOT support deletion, but Pilon = PASS (EBR = 0) - Illumina has FALSE confidence in a region which is NOT supported by the PacBio assembly alignment
            # Requirements: PB-MM2 DP = 0, Pilon-Tag == PASS, Within_PB_Alignment == False

        elif ( (PB_DP == 0) and (Ill_Filter == "PASS") and (PB_Supports_Del == False) ) :
            
            # Set EBR = np.nan (USED TO BE EBR = 0)
            listOfTuples.append( (ref_POS, np.nan, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["NAN_Counts"] += 1
            count_EBR_Outcomes_Dict["Outcome_8"] += 1 
 



        ### Outcomes 9 and 10 (Tandem Duplicated Regions in PB-MM2 Alignment)

        ## Outcome 9) PB-MM2 DP > 1, but Pilon = PASS (EBR = 0)
            # Requirements: PB-MM2 DP = 0, Pilon-Tag == PASS

        elif ( (PB_DP > 1) and (Ill_Filter == "PASS") ) : # 
            if (PB_ALT == Ill_ALT):
                # Set EBR = np.nan (USED TO BE EBR = 0)
                listOfTuples.append( (ref_POS, np.nan, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
                count_EBR_Outcomes_Dict["NAN_Counts"] += 1
                
                count_EBR_Outcomes_Dict["Outcome_9"] += 1  
                count_EBR_Outcomes_Dict["Outcome_9_DP2_NtAgree"] += 1

            else:
                # Set EBR = np.nan (USED TO BE EBR = 0)
                listOfTuples.append( (ref_POS, np.nan, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
                count_EBR_Outcomes_Dict["NAN_Counts"] += 1
                
                count_EBR_Outcomes_Dict["Outcome_9"] += 1  
                count_EBR_Outcomes_Dict["Outcome_9_DP2_NtDisagree"] += 1  

            
        ## Outcome 10) PB-MM2 DP > 1, but Pilon != PASS (EBR = 0)
            # Requirements: PB-MM2 DP = 0, Pilon-Tag != PASS

        elif ( (PB_DP > 1) and (Ill_Filter != "PASS") ) : # 

            # Set EBR = 0
            #listOfTuples.append( (ref_POS, 0, Ill_Filter))
            
            listOfTuples.append( (ref_POS, np.nan, Ill_Filter, Ill_TD, Ill_DP, Ill_MQ) )
            count_EBR_Outcomes_Dict["NAN_Counts"] += 1
            count_EBR_Outcomes_Dict["Outcome_10"] += 1             
    
        
        else: # Both technologies disagree BUT CANNOT BE CLASSIFIED
            print("Unknown Disagreement!!!! At position:", ref_POS)
            listOfTuples.append( (ref_POS, 0, Ill_Filter))
            
            count_EBR_Outcomes_Dict["Outcome_11_Unknown"] += 1
            
            
    #print("Breakdown of EBR Outcomes")
    #print(count_EBR_Outcomes_Dict)
    
    return listOfTuples, count_EBR_Outcomes_Dict




### Define function for output of a BED file from a NP array with values for each basepair position ###

# BED format specifications: https://useast.ensembl.org/info/website/upload/bed.html

def convert_GenomeNParray_To_BED_DF(input_GenomeNParray, genomeChrom = "NC_000962.3"):
    """ """
    last_Score = input_GenomeNParray[0]

    startOfRegion = 0
    listOfBED_Tuples = []
    RegionCounter = 1

    for RefPos_0based in (np.arange(len(input_GenomeNParray))):

        EBR_Score = input_GenomeNParray[RefPos_0based]

        if EBR_Score != last_Score:

            endOfRegion = RefPos_0based
            lengthOfRegion = endOfRegion - startOfRegion 

            BED_EntryTuple = (genomeChrom, startOfRegion, endOfRegion, f"Region{RegionCounter}_Length_{lengthOfRegion}_bp", last_Score,)
            listOfBED_Tuples.append(BED_EntryTuple)

            RegionCounter += 1

            #print(f"{H37rv_ChrName}, {startOfRegion}, {RefPos_0based}, {lengthOfRegion}_bp, {last_Score}, .")

            startOfRegion = RefPos_0based 

            #1 Output the last range
            #2 Store the new score    

        last_Score = EBR_Score #2 Store the new score   

        
        
    endOfRegion = RefPos_0based + 1
    lengthOfRegion = endOfRegion - startOfRegion 

    BED_EntryTuple = (genomeChrom, startOfRegion, endOfRegion, f"Region{RegionCounter}_Length_{lengthOfRegion}_bp", last_Score)
    listOfBED_Tuples.append(BED_EntryTuple)       

    BED_DF = pd.DataFrame(listOfBED_Tuples)
    
    BED_DF.columns = ["chrom", "chromStart", "chromEnd", "name", "score" ]
    
    
    return BED_DF


def convert_EBR_Array_To_TSV(input_EBR_Array, EBR_Output_TSV_PATH):
    
    input_EBR_TSV_DF = pd.DataFrame(input_EBR_Array)

    input_EBR_TSV_DF.columns = ["EmpiricalBasePairRecall"]
    input_EBR_TSV_DF["H37rv_RefPos_0based"] = input_EBR_TSV_DF.index
    input_EBR_TSV_DF["H37rv_RefPos_1based"] = input_EBR_TSV_DF["H37rv_RefPos_0based"] + 1

    input_EBR_TSV_DF = input_EBR_TSV_DF[["H37rv_RefPos_0based", "H37rv_RefPos_1based", "EmpiricalBasePairRecall"]  ]     

    input_EBR_TSV_DF.to_csv(EBR_Output_TSV_PATH, sep = "\t", index = False)
    
    print("TSV output to:", EBR_Output_TSV_PATH)
    
    return input_EBR_TSV_DF







print('Number of arguments:', len(sys.argv) -1, 'arguments.')
print( 'Argument List:', str(sys.argv) )
print("--------------------------------------------------------------\n")


def main():


    i_MM2_AtoRef_PP_BAM_mpileup_call_out_PATH = sys.argv[1]
    i_NucDiff_SVs_Filtered_DeletedPositionsOnly_GFF_PATH = sys.argv[2]
    i_Ill_Pilon_VCF_PATH = sys.argv[3]
    o_EBR_Agreement_DF_TSV_PATH = sys.argv[4]
    o_EBR_Outcome_Breakdown_Dict_JSON = sys.argv[5]
    o_EBR_BED_DF_PATH = sys.argv[6]
    o_EBR_BED_AmbOnly_DF_PATH = sys.argv[7]


    #print(i_MM2_AtoRef_PP_BAM_mpileup_call_out_PATH)


    print("STEP 1: Parsing PacBio Mpileup VCF")
    i_PB_mpileup_Dict = parse_PB_mpileup_ToDict(i_MM2_AtoRef_PP_BAM_mpileup_call_out_PATH)
    print("STEP 1: Finished Parsing PacBio Mpileup VCF")
    print("--------------------------------------------------------------\n")


    print("Step 2: Parsing PacBio NucDiff SVs GFF")
    i_PB_NucDiff_DeletionsDict = parse_NucDiff_SVs_DeletedPositions_FromGFF(i_NucDiff_SVs_Filtered_DeletedPositionsOnly_GFF_PATH)
    print("Step 2: Finished Parsing PacBio NucDiff SVs GFF")
    print("--------------------------------------------------------------\n")


    print("Step 3: Parsing Illumina Pilon VCF")
    i_Ill_Pilon_VCF_Dict = parse_Ill_PilonVCF_ToDict(i_Ill_Pilon_VCF_PATH)
    print("Step 3: Finished Parsing Illumina Pilon VCF")
    print("--------------------------------------------------------------\n")



    print("Step 4: Calculating Empirical Base Pair Recall (EBR) for all reference positions")
    print("--------------------------------------------------------------\n")
    i_ListOfTuples, i_count_EBR_Outcomes_Dict  = get_AgreementTuple_Between_PB_and_Ill_V7(i_PB_mpileup_Dict, i_Ill_Pilon_VCF_Dict, i_PB_NucDiff_DeletionsDict)

    print("Step 4: Finished calculating Empirical Base Pair Recall (EBR) for all reference positions")


    i_PBvsIll_Agreement_DF = pd.DataFrame(i_ListOfTuples)
    i_PBvsIll_Agreement_DF.columns = ["POS", "Agreement", "Ill_Pilon_Tag", "Ill_Pilon_TD", "Ill_Pilon_DP", "Ill_Pilon_MQ"]
    i_PBvsIll_Agreement_DF = i_PBvsIll_Agreement_DF[["Agreement", "Ill_Pilon_Tag", "Ill_Pilon_TD", "Ill_Pilon_DP", "Ill_Pilon_MQ"]]
    

    print("Outputting EBR DF to TSV format:")
    i_PBvsIll_Agreement_DF.to_csv(o_EBR_Agreement_DF_TSV_PATH, sep="\t", index = False)
    
    print("Outputting EBR Outcomes Dict to JSON format:")    
    with open(o_EBR_Outcome_Breakdown_Dict_JSON, "w") as outfile: 
    	json.dump(i_count_EBR_Outcomes_Dict, outfile) 



    
    i_EBR_Array = i_PBvsIll_Agreement_DF["Agreement"].fillna("Ambiguous").values

    i_EBR_BED_DF = convert_GenomeNParray_To_BED_DF(i_EBR_Array)

    i_EBR_BED_DF_AMB_ONLY = i_EBR_BED_DF[ i_EBR_BED_DF["score"] == "Ambiguous" ]


    i_EBR_BED_DF.to_csv(o_EBR_BED_DF_PATH,
                           sep = "\t",
                           index = False,
                           header = False)

    i_EBR_BED_DF_AMB_ONLY.to_csv(o_EBR_BED_AmbOnly_DF_PATH,
                           sep = "\t",
                           index = False,
                           header = False)



main()


