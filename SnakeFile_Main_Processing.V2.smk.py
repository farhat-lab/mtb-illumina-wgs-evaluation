# Snakemake - Run all core steps for Illumina and PacBio sequencing analysis
### Maximillian Marin (mgmarin@g.harvard.edu)



### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]

GATK4_PATH = config["GATK4_PATH"]


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

        
    if SampleID_i in input_SampleIDs_WithPacBio: # If PacBio WGS is provided

        PacBio_FQs_PATH = row["PacBio_FQ_PATH"]
        PacBio_FastQ_Files_List = PacBio_FQs_PATH.split(";")

        SampleID_To_PacBio_FQs_Dict[SampleID_i] = PacBio_FastQ_Files_List




# Define lists of aligners and variant callers to benchmark
listOf_Aligners = ["bwa-mem", "minimap2", "bowtie2"]  



listOf_VariantCallers = ["Pilon", "mpileup-call", "Varscan2"]

#listOf_VariantCallers = ["Pilon", "mpileup-call", "Varscan2", "HaplotypeCaller.2DCNN_scored", "HaplotypeCaller.HardFilt"]

listOf_VariantCallers = ["Pilon", "mpileup-call", "Varscan2", "HaplotypeCaller.2DCNN_scored", "HaplotypeCaller.HardFilt" ]

#listOf_VariantCallers = ["Pilon", "mpileup-call", "HaplotypeCaller.2DCNN_scored", "HaplotypeCaller.HardFilt" ]





#listOf_Aligners = ["bowtie2"]  

#listOf_VariantCallers = ["HaplotypeCaller.2DCNN_scored"]


# Define what metric to filter variants on for generating Precision-Recall curves of variant calling 

HapPy_VarCallers_FilterOn_QUAL = ["mpileup-call", "Varscan2"] #"HaplotypeCaller.HardFilt",
HapPy_VarCallers_FilterOn_MQ = ["Pilon"]
HapPy_VarCallers_FilterOn_CNN2D = ["HaplotypeCaller.2DCNN_scored"]
HapPy_VarCallers_FilterOn_QD = ["HaplotypeCaller.HardFilt"]



listOf_MaskSchemes_ToTest = ["NoRegionsRemoved", "PLCRegionsRemoved", "Pmap_K50E4_Below1Removed", "NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp", ]



##### Import Snakemake rules from SMK files #####
#include: "Snakemake_Rules_And_Config/1_PacBio_Mtb_AssemblyAndAnalysis.smk"
#include: "Snakemake_Rules_And_Config/2_Illumina_Mtb_Mtb_Analysis.smk"
#include: "Snakemake_Rules_And_Config/3_PacBio_And_Illumina_CombinedSteps.smk"
#include: "Snakemake_Rules_And_Config/4_SequencingBiasAnalysis.smk"
#include: "Snakemake_Rules_And_Config/5_Preprocessing_PB_Alignments_For_PhylogenyBuilding.smk"
include: "Snakemake_Rules_And_Config/7_Benchmark.All.Aligner.VariantCaller.Combos.smk"

# include: "Snakemake_Rules_And_Config/6_Evalulate_Masking_RLC_Regions.smk"

rule all:
    input:
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam",
               sampleID_WiIll = input_SampleIDs_WithIllumina,
               aligner_ID = listOf_Aligners),

        expand(output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.FilterTagged.vcf",
               sampleID_WiIll = input_SampleIDs_WithIllumina,
               aligner_ID = listOf_Aligners,
               variantCaller_ID = listOf_VariantCallers),


## Aligner-VariantCaller benchmarking ##

        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_{MaskScheme_ID}/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ.summary.csv",
               sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill,
               aligner_ID = listOf_Aligners,
               variantCaller_ID = HapPy_VarCallers_FilterOn_MQ,
               MaskScheme_ID = listOf_MaskSchemes_ToTest),


        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_{MaskScheme_ID}/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ.summary.csv",
               sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill,
               aligner_ID = listOf_Aligners,
               variantCaller_ID = HapPy_VarCallers_FilterOn_QUAL,
               MaskScheme_ID = listOf_MaskSchemes_ToTest),



        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_{MaskScheme_ID}/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ.summary.csv",
               sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill,
               aligner_ID = listOf_Aligners,
               variantCaller_ID = HapPy_VarCallers_FilterOn_CNN2D,
               MaskScheme_ID = listOf_MaskSchemes_ToTest),


        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_{MaskScheme_ID}/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ.summary.csv",
               sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill,
               aligner_ID = listOf_Aligners,
               variantCaller_ID = HapPy_VarCallers_FilterOn_QD,
               MaskScheme_ID = listOf_MaskSchemes_ToTest),











