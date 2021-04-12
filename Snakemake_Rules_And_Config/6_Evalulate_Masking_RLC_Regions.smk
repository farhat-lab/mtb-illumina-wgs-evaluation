

# Perform filtering of Pilon variant calls (1 to 15 bp Variants)

rule RemoveRefinedLowConfidenceRegions_From_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCF_PassOnly_Lengths_1to15bp:
    input:
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        RLC_Regions_BED = "references/210219.RefinedLCRs.bed"
    output:
        pilon_VCF_SNPsINDELs_1to15bp_RefinedLowConfidenceRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.RLCRegionsRemoved.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bedtools intersect -header -v -a {input.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -b {input.RLC_Regions_BED} -wa > {output.pilon_VCF_SNPsINDELs_1to15bp_RefinedLowConfidenceRegionsRemoved} "





rule Happy_AmbRemoved_VCeval_T_PB_G3PP_MM2_Paftools_Vs_Q_Ill_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCs_SNPsINDELs_Lengths_1to15bp_RefinedLowConfidenceRegionsRemoved:
    input:
        H37rv_FA = refGenome_FA_PATH,
        pilon_VCF_SNPsINDELs_1to15bp_RefinedLowConfidenceRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.RLCRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_RefinedLowConfidenceRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.pilon_VCF_SNPsINDELs_1to15bp_RefinedLowConfidenceRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{wildcards.sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_RefinedLowConfidenceRegionsRemoved/Hap.py.{wildcards.sampleID_Wi_Ill_And_PB} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "


###################################################################################################
############################################## END ################################################
###################################################################################################

