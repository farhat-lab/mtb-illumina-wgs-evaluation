

###########################################################
######### Phylogeny Building Pre-processing Steps #########
###########################################################


### Merging VCFs w/ ALT and REF alleles represented across all samples###


## 1) Get list of SNP positions from each sample VCF (PASS ONLY, NO MQ FILTERING (MQ >=0))

### cut -f 1,2 <__>

## 2) Merge w/ cat, sort, uniq to create a deduplicated list of all SNP positions

## 3) Convert and compress VCF to bcf.gz (and index)

## 4) 



##### Merge SNPs for Phylogeny generation (Illumina WGS) #####

rule getAll_SNPpositions_Pilon_SNPsINDELsOnly_VarCalling_VCF_ToH37rv:
    input:
        pilon_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf",
    output:
        pilon_VCF_AllSNPpositions_TSV = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.AllSNPpositions.PilonPASS.SNPs.tsv",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view --types snps -f PASS {input.pilon_VCF} | cut -f 1,2 | grep -v '#' > {output.pilon_VCF_AllSNPpositions_TSV} \n"



awk_AddH37rvChromoColumn_String = '{print "NC_000962.3\t"$0}'

rule combineAll_SNPpositions_Pilon_SNPsINDELsOnly_VarCalling_VCF_ToH37rv:
    input:
        expand(output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.AllSNPpositions.PilonPASS.SNPs.tsv", sampleID_WiIll=input_SampleIDs_WithIllumina),
    output:
        AllSample_AllSNPpositions_TSV = output_Dir + "/MergingSNPS_AcrossAllSamples/AllSNPpositions.PilonPASS.SNPs.Union.AllSamples.tsv"
    shell:
        "cat {input} | sort -k 2n | uniq > {output}"



rule Filter_Pilon_VarCalling_To_OnlySNPpositionsInUnionOfAllSNPs:
    input:
        pilon_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf",
        AllSample_AllSNPpositions_TSV = output_Dir + "/MergingSNPS_AcrossAllSamples/AllSNPpositions.PilonPASS.SNPs.Union.AllSamples.tsv",
    output:
        pilon_BCF_GZ = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.bcf.gz",
        pilon_BCF_SNPsInAllSamples = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.OnlySNPpositionsIn.UnionOfAllSamples.bcf.gz",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view {input.pilon_VCF} -O b -o {output.pilon_BCF_GZ} \n"
        "bcftools index {output.pilon_BCF_GZ} \n"
        "bcftools view {output.pilon_BCF_GZ} -f PASS "
        " -R {input.AllSample_AllSNPpositions_TSV} -O b -o {output.pilon_BCF_SNPsInAllSamples} \n"
        "bcftools index {output.pilon_BCF_SNPsInAllSamples} \n"

#######################################################################







##### Merge SNPs for Phylogeny generation (PacBio) #####

rule getAll_SNPpositions_PBMM2_mpileup_VarCalling_VCF_ToH37rv:
    input:
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.KeepAllPositions.vcf"
    output:
        PBMM2_mpileup_VCF_AllSNPpositions_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.AllSNPpositions.tsv",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view --types snps {input.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF} | cut -f 1,2 | grep -v '#' > {output.PBMM2_mpileup_VCF_AllSNPpositions_TSV} \n"



rule combineAll_SNPpositions_PBMM2_mpileup_VarCalling_VCF_ToH37rv:
    input:
        expand(output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.AllSNPpositions.tsv", sampleID_Wi_Ill_And_PB=input_SampleIDs_With_PB_And_Ill),           
    output:
        AllSample_AllSNPpositions_PBMM2_mpileup_TSV = output_Dir + "/MergingSNPS_AcrossAllSamples/AllSNPpositions.PBMM2.mpileup.SNPs.Union.AllSamples.tsv"
    shell:
        "cat {input} | sort -k 2n | uniq > {output}"



rule Filter_PBMM2_mpileup_VarCalling_To_OnlySNPpositionsInUnionOfAllSNPs:
    input:
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
        AllSample_AllSNPpositions_PBMM2_mpileup_TSV = output_Dir+ "/MergingSNPS_AcrossAllSamples/AllSNPpositions.PBMM2.mpileup.SNPs.Union.AllSamples.tsv",
    output:
        MM2_Isolate_To_Ref_AllPositions_BCF_GZ = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.AllPositions.bcf.gz",
        MM2_Isolate_To_Ref_SNPs_BCF_GZ_SNPsInAllSamples = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.mpileup.call.SNPs.Union.AllSamples.bcf.gz",

    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view {input.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF} -O b -o {output.MM2_Isolate_To_Ref_AllPositions_BCF_GZ} \n"
        "bcftools index {output.MM2_Isolate_To_Ref_AllPositions_BCF_GZ} \n"
        "bcftools view {output.MM2_Isolate_To_Ref_AllPositions_BCF_GZ} "
        " -R {input.AllSample_AllSNPpositions_PBMM2_mpileup_TSV} -O b -o {output.MM2_Isolate_To_Ref_SNPs_BCF_GZ_SNPsInAllSamples} \n"
        "bcftools index {output.MM2_Isolate_To_Ref_SNPs_BCF_GZ_SNPsInAllSamples} \n"

#######################################################################





