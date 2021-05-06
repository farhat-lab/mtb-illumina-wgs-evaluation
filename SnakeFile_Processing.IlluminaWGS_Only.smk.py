# Snakemake - Run all core steps for Illumina and PacBio sequencing analysis
### Maximillian Marin (mgmarin@g.harvard.edu)



### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]


# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]

# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')

input_DataInfo_DF_With_PacBio_WGS = input_DataInfo_DF[  input_DataInfo_DF["PacBio_FQ_PATH"] != "None" ]
input_DataInfo_DF_With_Illumina_WGS = input_DataInfo_DF[ input_DataInfo_DF["Illumina_PE_FQs_PATH"] != "None" ]
input_DataInfo_DF_With_Both_Illumina_And_PacBio_WGS = input_DataInfo_DF[ (input_DataInfo_DF["Illumina_PE_FQs_PATH"] != "None") & (input_DataInfo_DF["PacBio_FQ_PATH"] != "None") ]


# Create a python list of Sample IDs

input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )
input_SampleIDs_WithPacBio = list( input_DataInfo_DF_With_PacBio_WGS["SampleID"].values )
input_SampleIDs_WithIllumina = list( input_DataInfo_DF_With_Illumina_WGS["SampleID"].values )
input_SampleIDs_With_PB_And_Ill = list( input_DataInfo_DF_With_Both_Illumina_And_PacBio_WGS["SampleID"].values )


SampleIDTo_Illumina_PE_FQ1_Dict = {}
SampleIDTo_Illumina_PE_FQ2_Dict = {}

SampleID_To_PacBio_FQs_Dict = {}


# Iterate over each sample's row and define path to input FQ files for each sequencing technology

for idx, row in input_DataInfo_DF.iterrows():
    
    SampleID_i = row["SampleID"]
    
    if SampleID_i in input_SampleIDs_WithIllumina: # If Illumina WGS is provided

        Illumina_FQs_PATH = row["Illumina_PE_FQs_PATH"]
        Ill_FastQ_Files_List = Illumina_FQs_PATH.split(";")
        

        if len(Ill_FastQ_Files_List) == 2:

            FQ_1_PATH, FQ_2_PATH = Ill_FastQ_Files_List
            SampleIDTo_Illumina_PE_FQ1_Dict[SampleID_i] = FQ_1_PATH
            SampleIDTo_Illumina_PE_FQ2_Dict[SampleID_i] = FQ_2_PATH

        
#    if SampleID_i in input_SampleIDs_WithPacBio: # If PacBio WGS is provided
#
#        PacBio_FQs_PATH = row["PacBio_FQ_PATH"]
#        PacBio_FastQ_Files_List = PacBio_FQs_PATH.split(";")

#        SampleID_To_PacBio_FQs_Dict[SampleID_i] = PacBio_FastQ_Files_List







##### Import Snakemake rules from SMK files #####
#include: "Snakemake_Rules_And_Config/1_PacBio_Mtb_AssemblyAndAnalysis.smk"
include: "Snakemake_Rules_And_Config/2_Illumina_Mtb_Mtb_Analysis.smk"
include: "Snakemake_Rules_And_Config/3_PacBio_And_Illumina_CombinedSteps.smk"
#include: "Snakemake_Rules_And_Config/4_SequencingBiasAnalysis.smk"
#include: "Snakemake_Rules_And_Config/5_Preprocessing_PB_Alignments_For_PhylogenyBuilding.smk"
include: "Snakemake_Rules_And_Config/6_Evalulate_Masking_RLC_Regions.smk"

rule all:
    input:
        output_Dir + "/multiqc_Reports/multiqc_AllDirs.html",
        output_Dir + "multiqc_Reports/multiqc_Trimmomatic.html",
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiIll}.IllPE.H37rv.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.PassOnly.SNPsOnly.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),

        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/bcftools_mpileup_call_IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.bcftools.call.IllPE.H37rv.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),

        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.depth.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.depth.averaged.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),


        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir +"/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.RLCRegionsRemoved.vcf", sampleID_WiIll = input_SampleIDs_WithIllumina),






        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/{sampleID_WiIll}_F2.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_WiIll}/LineageCalling/LineageCall_IlluminaWGS_AlignTo_H37rv/{sampleID_WiIll}.IlluminaWGS.Pilon.lineage_call.Coll.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),



