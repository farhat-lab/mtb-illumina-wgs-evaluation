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

        
    if SampleID_i in input_SampleIDs_WithPacBio: # If PacBio WGS is provided

        PacBio_FQs_PATH = row["PacBio_FQ_PATH"]
        PacBio_FastQ_Files_List = PacBio_FQs_PATH.split(";")

        SampleID_To_PacBio_FQs_Dict[SampleID_i] = PacBio_FastQ_Files_List







##### Import Snakemake rules from SMK files #####
include: "Snakemake_Rules_And_Config/1_PacBio_Mtb_AssemblyAndAnalysis.smk"
include: "Snakemake_Rules_And_Config/2_Illumina_Mtb_Mtb_Analysis.smk"
include: "Snakemake_Rules_And_Config/3_PacBio_And_Illumina_CombinedSteps.smk"
include: "Snakemake_Rules_And_Config/4_SequencingBiasAnalysis.smk"
include: "Snakemake_Rules_And_Config/5_Preprocessing_PB_Alignments_For_PhylogenyBuilding.smk"
# include: "Snakemake_Rules_And_Config/6_Evalulate_Masking_RLC_Regions.smk"

rule all:
    input:
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/Nanoplot_QC/{sampleID_WiPacBio}.subreads.NanoPlot/NanoPlot-report.html", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.ReadLengths.tsv", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.depth.txt", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.depth.averaged.txt", sampleID_WiPacBio = input_SampleIDs_WithPacBio),


        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta", sampleID_WiPacBio = input_SampleIDs_WithPacBio),


        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/{sampleID_WiIll}_F2.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/F2_Calculation/{sampleID_WiPacBio}_F2.txt", sampleID_WiPacBio = input_SampleIDs_WithPacBio),



        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/FastANI_OutputDirs/FlyeAssembly_I3_PBonly_FastANI/{sampleID_WiPacBio}.flyeassembly.I3.AssemblyToH37rv.FastANI.txt", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_WiPacBio}/Prokka/Flye_I3_PB_DeNovo_Prokka_Anno/{sampleID_WiPacBio}.fna", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/FastANI_OutputDirs/Flye_I3_PilonPolished_FastANI/{sampleID_Wi_Ill_And_PB}_I3_PP_AssemblyToH37rv.FastANI.txt", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Prokka/Flye_I3_PB_PilonPolished_Prokka_Anno/{sampleID_Wi_Ill_And_PB}.fna", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        


        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_CoscollaRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),                     
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_SNPsINDELs_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp_V2/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),   

        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.mpileup.call.SNPs.Union.AllSamples.bcf.gz", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.OnlySNPpositionsIn.UnionOfAllSamples.bcf.gz", sampleID_WiIll = input_SampleIDs_WithIllumina),


        expand(output_Dir + "/{sampleID_WiIll}/LineageCalling/LineageCall_IlluminaWGS_AlignTo_H37rv/{sampleID_WiIll}.IlluminaWGS.Pilon.lineage_call.Coll.txt", sampleID_WiIll = input_SampleIDs_WithIllumina),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/LineageCalling/LineageCall_Flye_I3_PP_MM2_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.AtoRef.FlyeI3_PP.lineage_call.Coll.txt", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_WiPacBio}/LineageCalling/LineageCall_Flye_I3_PBonly_MM2_AlignTo_H37rv/{sampleID_WiPacBio}.AtoRef.Flye_I3_PBonly.lineage_call.Coll.txt", sampleID_WiPacBio = input_SampleIDs_WithPacBio),



        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam.depth.txt", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.depth.txt", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),


## RLC Extra Step
        #expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_RefinedLowConfidenceRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv", sampleID_Wi_Ill_And_PB = input_SampleIDs_With_PB_And_Ill),  


