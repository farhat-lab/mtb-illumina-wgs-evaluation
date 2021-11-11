###### Illumina Variant Calling & Benchmarking - Including multiple Aligner-Caller combinations  #####


### Adding multiple Aligner-Caller combinations  #####




# For each combination of a) Aligner and b) Variant Caller, we want to do the following

### 1) Align with choosen aligner (BWA-mem, Minimap2, URMAP, and bowtie2)

### 2) Use PICARD to remove duplicates

### 3) Run variant calling with choosen software (Bcftools-mpileup, VarScan2, GATK-HaplotypeCaller)




# Proposed Directory Organization

### 
### AlignerAndCaller_Benchmarking/
###      [Aligner_A]/
###               Alignments_[Aligner_A]/
###               VarCall_[VariantCaller_A]/
###               VarCall_[VariantCaller_B]/
###               Happy_Benchmarking_[Aligner_A]_[VariantCaller_A]/
###               Happy_Benchmarking_[Aligner_A]_[VariantCaller_B]/
###
###      [Aligner_B]/
###               Alignments_[Aligner_B]/
###               VarCall_[VariantCaller_A]/
###               VarCall_[VariantCaller_B]/
###               Happy_Benchmarking_[Aligner_B]_[VariantCaller_A]/
###               Happy_Benchmarking_[Aligner_B]_[VariantCaller_B]/




### 
### AlignerAndCaller_Benchmarking/
### bwa-mem/
###         Alignments/
###         VarCall_Pilon/
###         VarCall_bcftools-mpileup/
###         VarCall_Varscan2/
###         VarCall_GATK-HaplotypeCaller-HardFilt/
###         Happy_Benchmarking_bwa-mem_Pilon/
###         Happy_Benchmarking_bwa-mem_bcftools-mpileup/
###         Happy_Benchmarking_bwa-mem_Varscan2/
###         Happy_Benchmarking_bwa-mem_GATK-HaplotypeCaller-HardFilt/
###


#### Alignment of Illumina Paired End Reads ###


# 1) bwa-mem (Default parameters)

rule bwa_map_IllPE_AlignTo_H37rv_AllAligners:
    input:
        fa = refGenome_FA_PATH,
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/bwa-mem/Alignments_bwa-mem/{sampleID_WiIll}.IllPE.bwa-mem.sam"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_WiIll}\tSM:{sampleID_WiIll}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.fa} {input.fq1_trimmed} {input.fq2_trimmed} > {output}"


# 2) minimap2 (Using default parameters for short reads)


# 3) bowtie2 (Using default parameters for short reads)


# 4) URMAP (Using default parameters for short reads)








### Alignment Processing - Reformating and removal of duplicates steps ###

##### NOTE: These will be the same across ALL ALIGNERS! #####


rule samtools_ViewAndSort_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.sam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} -m 4G | samtools sort -m 4G - > {output}"



rule samtools_index_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam.bai"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"


#####################################
#### PICARD (remove duplicates) #####
#####################################

rule picard_RemoveDup_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    output:
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_BwaMEM_Duprem_METRICS = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.metrics",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx6g MarkDuplicates I={input} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"


rule samtools_index_IllPE_Duprem_AlignTo_H37rv:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",    output:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.bai",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"






# Using Pilon with "--fix all,breaks --minmq 1 --mindepth 5" 

rule pilon_VarCalling_IllPE_AlignTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks:
    input:
        Ref_fa = refGenome_FA_PATH,
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_BwaMEM_Duprem_BAM_BAI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/bwa-mem/Alignments_bwa-mem/{sampleID_WiIll}.IllPE.bwa-mem.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/bwa-mem/VarCall_Pilon/{sampleID_WiIll}.IllPE.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx12g --genome {input.Ref_fa} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_WiIll}.IllPE.H37rv"
        " --outdir {params.Pilon_OutputDir_PATH} --fix all,breaks --vcf --changes --tracks --minmq 1 --mindepth 5"







### Filtering and masking VCFs - for ALL variant callers ###


### Filter out all INDELs with length greater than 15 bp 
rule Filter_Pilon_VCF_minMQ_1_minDP_5_Fix_All_Breaks_RemoveIndelsGreaterThan15bp:
    input:
        pilon_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.vcf",
    output:
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view --types snps,indels -i 'abs(strlen(ALT)-strlen(REF))<=15' -f PASS {input.pilon_VCF} > {output.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only}"


# Perform masking of PLC and Pmap-K50E4 < 1 (Pileup Mappability)

rule RemoveCoscollaRegions_From_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCF_PassOnly_Lengths_1to15bp:
    input:
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        Coscolla_etal_Regions_BED = "References/Mtb_H37Rv_MaskingSchemes/201027_Mtb_H37rv_pLC_Regions_CoscollaExcludedGenes.bed"
    output:
        pilon_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bedtools intersect -header -v -a {input.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -b {input.Coscolla_etal_Regions_BED} -wa > {output.pilon_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} "



rule Remove_Pmap_K50E4_Below1_From_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCF_PassOnly_Lengths_1to15bp:
    input:
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        Pmap_K50E4_Below1_Regions_BED = "References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed",
    output:
        pilon_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bedtools intersect -header -v -a {input.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -b {input.Pmap_K50E4_Below1_Regions_BED} -wa > {output.pilon_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} "






##### Running Hap.py Benchmarking - For ALL Aligner-VariantCaller combinations ######

##### NOTE: These rules will be the same across ALL ALIGNER-VariantCaller combinations! #####

rule Happy_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved:
    input:
        H37rv_FA = refGenome_FA_PATH,
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{wildcards.sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_NoRegionsRemoved/Hap.py.{wildcards.sampleID_Wi_Ill_And_PB} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "




rule Happy_SNPsINDELs_Lengths_1to15bp_CoscollaRegionsRemoved:
    input:
        H37rv_FA = refGenome_FA_PATH,
        pilon_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_CoscollaRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.pilon_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{wildcards.sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_CoscollaRegionsRemoved/Hap.py.{wildcards.sampleID_Wi_Ill_And_PB} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "





rule Happy_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed:
    input:
        H37rv_FA = refGenome_FA_PATH,
        pilon_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.pilon_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} -r {input.H37rv_FA} " 
        "-o {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{wildcards.sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_Pmap_K50E4_Below1Removed/Hap.py.{wildcards.sampleID_Wi_Ill_And_PB} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "



rule Happy_All_WiStratificationBy_Pmap_and_SV_4sets_V2:
    input:
        H37rv_FA = refGenome_FA_PATH,
        pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_SNPsINDELs_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp_V2/Hap.py.{sampleID_Wi_Ill_And_PB}.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.pilon_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/PBMM2_Paftools_GroundTruthVCF_Evaluations_V3_minMQ_1_minDP_5_Fix_All_Breaks_AmbRegionsRemoved/{wildcards.sampleID_Wi_Ill_And_PB}_Happy_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_VCs_SNPsINDELs_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp_V2/Hap.py.{wildcards.sampleID_Wi_Ill_And_PB} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} "
        " --stratification {input.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} --roc-regions HighPmap_NoSV,LowPmap_NoSV,LowPmap_WiSV,HighPmap_WiSV -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "














