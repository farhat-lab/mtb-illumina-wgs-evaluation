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



# Example of directory structure for bwa-mem:
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

rule bwa_map_IllPE:
    input:
        H37rv_FA = refGenome_FA_PATH,
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/bwa-mem/Alignments_bwa-mem/{sampleID_WiIll}.IllPE.bwa-mem.sam"
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_WiIll}\tSM:{sampleID_WiIll}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.H37rv_FA} {input.fq1_trimmed} {input.fq2_trimmed} > {output}"


# 2) minimap2 (Using default parameters for short reads)

rule MM2_SR_IllPE:
    input:
        H37rv_FA = refGenome_FA_PATH,
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/minimap2/Alignments_minimap2/{sampleID_WiIll}.IllPE.minimap2.sam"
    threads: 8
    params:
        rg=r"@RG\tID:{sampleID_WiIll}\tSM:{sampleID_WiIll}"
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "minimap2 -R '{params.rg}' -ax sr -t {threads} {input.H37rv_FA} {input.fq1_trimmed} {input.fq2_trimmed} > {output} "

  
# 3) bowtie2 (Using default parameters for short reads)

# Step 1: Create index of reference
rule make_Ref_IDX_bowtie2:
    input:
        H37rv_FA = refGenome_FA_PATH,
    output:
        H37rv_IDX_BT2 = ".".join(refGenome_FA_PATH.split(".")[:-1]) + ".1.bt2"
    params:
        BT2_IDX_Out_Prefix = ".".join(refGenome_FA_PATH.split(".")[:-1])
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bowtie2-build {input} {params.BT2_IDX_Out_Prefix}"

# Step 2: Align reads to reference

rule bowtie2_IllPE:
    input:
        H37rv_IDX_BT2 = ".".join(refGenome_FA_PATH.split(".")[:-1]) + ".1.bt2",
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/bowtie2/Alignments_bowtie2/{sampleID_WiIll}.IllPE.bowtie2.sam"
    threads: 8
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    params:
        BT2_IDX_Prefix = ".".join(refGenome_FA_PATH.split(".")[:-1]) ,
        rg_SM="SM:{sampleID_WiIll}"
        #rg=r"@RG\tID:{sampleID_WiIll}\tSM:{sampleID_WiIll}"
    shell:
        "bowtie2 --rg-id {wildcards.sampleID_WiIll} --rg '{params.rg_SM}' -p {threads} -x {params.BT2_IDX_Prefix} -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} > {output} "



# 4) URMAP (Using default parameters for short reads)

# Step 1: Create index of reference
rule make_Ref_IDX_URMAP:
    input:
        H37rv_FA = refGenome_FA_PATH,
    output:
        H37rv_ID_UFI = output_Dir + "/URMAP_Ref_IDX/InputRef.H37Rv.urmap.ufi"
    shell:
        "urmap -make_ufi {input} -output {output}"

# Step 2: Align reads to reference
rule URMAP_IllPE:
    input:
        H37rv_ID_UFI = output_Dir + "/URMAP_Ref_IDX/InputRef.H37Rv.urmap.ufi",
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/urmap/Alignments_urmap/{sampleID_WiIll}.IllPE.urmap.sam"
    threads: 8
    shell:
        "urmap --threads {threads} -map2 {input.fq1_trimmed} -reverse {input.fq2_trimmed} -ufi {input.H37rv_ID_UFI} -samout {output}"





### Alignment Processing - Reformating and removal of duplicates steps ###

##### NOTE: These will be the same across ALL ALIGNERS! #####


rule samtools_ViewAndSort_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.sam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "samtools view -bS {input} -m 4G | samtools sort -m 4G - > {output}"



rule samtools_index_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam.bai"
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
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
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "picard -Xmx6g MarkDuplicates I={input} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"


rule samtools_index_IllPE_Deduped_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.bai",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "samtools index {input}"



rule samtools_rmdup_IllPE_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.bam"
    output:
        IllPE_BwaMEM_st_rmdp_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.rmdup.bam",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "samtools rmdup {input} {output} "


rule samtools_index_Of_IllPE_samtools_rmdup_AllAligners:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.rmdup.bam",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.rmdup.bam.bai",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "samtools index {input}"









###### Running Variant Callers ######




####  Pilon - Variant Calling  ####

# Pilon - with "--fix all,breaks --minmq 1 --mindepth 5" 

rule pilon_VariantCalling_IllPE:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_Duprem_BAM_BAI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Pilon/{sampleID_WiIll}.IllPE.{aligner_ID}.Pilon.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Pilon/{sampleID_WiIll}.IllPE.{aligner_ID}.Pilon.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Pilon/",
        Pilon_Output_Prefix = "{sampleID_WiIll}.IllPE.{aligner_ID}.Pilon"
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "pilon -Xmx12g --genome {input.H37rv_FA} --bam {input.IllPE_Duprem_BAM} --output {params.Pilon_Output_Prefix} "
        " --outdir {params.Pilon_OutputDir_PATH} --fix all,breaks --vcf --changes --tracks --minmq 1 --mindepth 5"


rule pilon_VariantCalling_IllPE_GetVars_FilterTaggedFile:
    input:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Pilon/{sampleID_WiIll}.IllPE.{aligner_ID}.Pilon.vcf",
    output:
        pilon_VCF_FilterTagged = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.Pilon.FilterTagged.vcf",   
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bcftools view --types snps,indels {input.pilon_VCF} > {output.pilon_VCF_FilterTagged}"



#### samtools mpileup ####

rule bcftools_Call_VariantCalling_IllPE:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_Duprem_BAM_BAI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.bai",
    output:
        mpileup_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_mpileup-call/{sampleID_WiIll}.IllPE.{aligner_ID}.mpileup-call.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bcftools mpileup -Ou -f {input.H37rv_FA} {input.IllPE_Duprem_BAM} "
        " | bcftools call -vm -O v -o {output.mpileup_VCF} "



# bcftools-call - Label "PASS" variants that do not pass certain HARD-Filtering thresholds

## NOTE: The resulting VCF can then be benchmarked for using HARD-filtering & for CNN based quality annotation

# Link to suggested standard filtering params: https://samtools.github.io/bcftools/howtos/variant-calling.html
# GATK Guide to variant filtration: https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2

rule bcftools_Call_VariantCalling_Filtration: # 
    input:
        H37rv_FA = refGenome_FA_PATH,
        mpileup_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_mpileup-call/{sampleID_WiIll}.IllPE.{aligner_ID}.mpileup-call.vcf",
    output:
        mpileup_FilterTagged_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_mpileup-call/{sampleID_WiIll}.IllPE.{aligner_ID}.mpileup-call.FilterTagged.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        'gatk VariantFiltration -R {input.H37rv_FA} -V {input.mpileup_VCF} -O {output.mpileup_FilterTagged_VCF} '
        ' --filter-name "RPB01"  --filter-expression "RPB < 0.1" '
        ' --filter-name "MQSB01"  --filter-expression "MQSB < 0.1" '
        ' --filter-name "DP10"  --filter-expression "DP < 10" '










#### VarScan2 (Using mpileup for initial processing)  ####


rule mpileup_VarScan2_VariantCalling_IllPE_SNPs_And_INDELs:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_Duprem_BAM_BAI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.bai",
    output:
        mpileup_VarScan2_SNPs_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.snps.vcf",
        mpileup_VarScan2_INDELs_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.indels.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "samtools mpileup -B -f {input.H37rv_FA} {input.IllPE_Duprem_BAM} "
        " | varscan mpileup2snp  --output-vcf 1  > {output.mpileup_VarScan2_SNPs_VCF} \n "
        ""
        "samtools mpileup -B -f {input.H37rv_FA} {input.IllPE_Duprem_BAM} "
        " | varscan mpileup2indel  --output-vcf 1  > {output.mpileup_VarScan2_INDELs_VCF} \n"




rule BGZIP_And_Tabix_VarScan2_SNPsONLY:
    input:
        mpileup_VarScan2_SNPs_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.snps.vcf",
    output:
        mpileup_VarScan2_SNPs_VCF_GZ = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.snps.vcf.gz",
        mpileup_VarScan2_SNPs_VCF_GZ_TBI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.snps.vcf.gz.tbi",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
         "bgzip --stdout {input.mpileup_VarScan2_SNPs_VCF} > {output.mpileup_VarScan2_SNPs_VCF_GZ}  \n"
         "tabix {output.mpileup_VarScan2_SNPs_VCF_GZ}"

rule BGZIP_And_Tabix_VarScan2_INDELsONLY:
    input:
        mpileup_VarScan2_INDELs_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.indels.vcf",
    output:
        mpileup_VarScan2_INDELs_VCF_GZ = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.indels.vcf.gz",
        mpileup_VarScan2_INDELs_VCF_GZ_TBI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.indels.vcf.gz.tbi",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
         "bgzip --stdout {input.mpileup_VarScan2_INDELs_VCF} > {output.mpileup_VarScan2_INDELs_VCF_GZ}  \n"
         "tabix {output.mpileup_VarScan2_INDELs_VCF_GZ}"


rule mergeVCFs_SNPandINDEL_VarScan2:
    input:
        mpileup_VarScan2_SNPs_VCF_GZ = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.snps.vcf.gz",
        mpileup_VarScan2_INDELs_VCF_GZ = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.indels.vcf.gz",
    output:
        mpileup_VarScan2_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        "bcftools concat -O v --allow-overlaps {input} > {output}"



rule VarScan2_Update_VCF_GQtoQUAL:
    input:
        mpileup_VarScan2_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.vcf",
    output:
        mpileup_VarScan2_GQtoQUAL_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.GQtoQUAL.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        "python Scripts/Varscan.VCF.Updater.GQtoQUAL.py --input {input} > {output}"



# VarScan2 - Label "PASS" variants that do not pass certain HARD-Filtering thresholds

## NOTE: The resulting VCF can then be benchmarked for using HARD-filtering & for CNN based quality annotation

# Link to suggested standard filtering params: https://samtools.github.io/bcftools/howtos/variant-calling.html
# GATK Guide to variant filtration: https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2

 # Note: This step is not done with GATK's VariantFiltration, instead it use's Varscan's filtration function
         #### This is b/c Varscan puts the variant annotations at the SAMPLE level, this makes it hard to use standard tools to filter the VCF.


rule VarScan2_VariantCalling_Filtration:
    input:
        mpileup_VarScan2_GQtoQUAL_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.GQtoQUAL.vcf",
    output:
        mpileup_VarScan2_FilterTagged_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_Varscan2/{sampleID_WiIll}.IllPE.{aligner_ID}.Varscan2.FilterTagged.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        "varscan filter {input} --min-coverage 10 --min-reads2 8 --min-avg-qual 20  --min-var-freq 0.8 --output-file {output}"










#### GATK - HaplotypeCaller ####

rule gatk_HaplotypeCaller_VariantCalling_IllPE:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        IllPE_Duprem_BAM_BAI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam.bai",
    output:
        gatk_HapCaller_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.UnfilteredOutput.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        "gatk HaplotypeCaller -R {input.H37rv_FA} -I {input.IllPE_Duprem_BAM} -O {output.gatk_HapCaller_VCF} --sample-ploidy 1"


rule gatk_CNNscore_2D:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_Duprem_BAM = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Alignments_{aligner_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.duprem.bam",
        gatk_HapCaller_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.UnfilteredOutput.vcf",
    output:
        gatk_HapCaller_CNN_2D_Scored_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.RawHaploidOut.vcf",
    #conda:
    #    "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    conda:
        "CondaEnvs/gatkcondaenv.4.1.7.0._MGMmod.yml"
    shell:
        "{GATK4_PATH} CNNScoreVariants -I {input.IllPE_Duprem_BAM} -V {input.gatk_HapCaller_VCF} "
        " -R {input.H37rv_FA} -O {output.gatk_HapCaller_CNN_2D_Scored_VCF} "
        " --tensor-type read_tensor --transfer-batch-size 8 --inference-batch-size 8 \n"



sed_ReplaceCommand_HapToDip = 's/GT:AD:DP:GQ:PL\t1/GT:AD:DP:GQ:PL\t1\/1/g'

rule convert_GT_HaploidToDiploid_gatk_CNNscore_2D_IllPE_AlignTo_H37rv: # Replace all GT of 1 to 1/1 in {SAMPLEID}.IllPE.GATK.2D_CNN_Filtered.vcf
    input:
        H37rv_FA = refGenome_FA_PATH,
        gatk_HapCaller_CNN_2D_Scored_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.RawHaploidOut.vcf",
    output:
        gatk_HapCaller_CNN_2D_Scored_DiploidGT_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        'sed "{sed_ReplaceCommand_HapToDip}" {input.gatk_HapCaller_CNN_2D_Scored_VCF} > {output.gatk_HapCaller_CNN_2D_Scored_DiploidGT_VCF} \n'
        

        #gatk_HapCaller_CNN_2D_Scored_DiploidGT_FilterTagged_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.CNN_2D_Filter.vcf"

        #'{GATK4_PATH} VariantFiltration -R {input.H37rv_FA} -V {output.gatk_HapCaller_CNN_2D_Scored_DiploidGT_VCF} -O {output.gatk_HapCaller_CNN_2D_Scored_DiploidGT_FilterTagged_VCF} '
        #'--filter-name "CNN_2D_Filter" --filter-expression "CNN_2D <= 0.0" '



# GATK-CNN - Label "FAIL" variants that do not pass certain HARD-Filtering thresholds

### OR Label "PASS" variants only

## NOTE: The resulting VCF can then be benchmarked for using HARD-filtering & for CNN based quality annotation


# VarScan2 - Label "PASS" variants that do not pass certain HARD-Filtering thresholds

## NOTE: The resulting VCF can then be benchmarked for using HARD-filtering & for CNN based quality annotation

# Link to suggested standard filtering params: https://samtools.github.io/bcftools/howtos/variant-calling.html
# GATK Guide to variant filtration: https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2

rule GATK_CNN_VariantCalling_Filtration: # 
    input:
        H37rv_FA = refGenome_FA_PATH,
        gatk_HapCaller_CNN_2D_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.vcf",
    output:
        gatk_HapCaller_CNN_2D_FilterTagged_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.FilterTagged.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        'gatk VariantFiltration -R {input.H37rv_FA} -V {input.gatk_HapCaller_CNN_2D_VCF} -O {output.gatk_HapCaller_CNN_2D_FilterTagged_VCF} '
        #'--filter-name "CNN_2D_Filter" --filter-expression "CNN_2D <= 0.0" '
        ' --filter-name "QD2"  --filter-expression "QD < 2.0"  '
        ' --filter-name "QUAL30"  --filter-expression "QUAL < 30.0" '
        ' --filter-name "SOR3"  --filter-expression "SOR > 3.0" '
        ' --filter-name "FS60"  --filter-expression "FS > 60.0" '
        ' --filter-name "DP10"  --filter-expression "DP < 10.0" '
        ' --filter-name "MQRankSum-12.5"  --filter-expression "MQRankSum < -12.5" '
        ' --filter-name "ReadPosRankSum-8"  --filter-expression "ReadPosRankSum < -8.0" '



rule GATK_StandardHardFilt_VariantCalling_Filtration: # 
    input:
        H37rv_FA = refGenome_FA_PATH,
        gatk_HapCaller_CNN_2D_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.2DCNN_scored/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.2DCNN_scored.vcf",
    output:
        gatk_HapCaller_HardFilt_FilterTagged_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_HaplotypeCaller.HardFilt/{sampleID_WiIll}.IllPE.{aligner_ID}.HaplotypeCaller.HardFilt.FilterTagged.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml" 
    shell:
        'gatk VariantFiltration -R {input.H37rv_FA} -V {input.gatk_HapCaller_CNN_2D_VCF} -O {output.gatk_HapCaller_HardFilt_FilterTagged_VCF} '
        ' --filter-name "QD2"  --filter-expression "QD < 2.0"  '
        ' --filter-name "QUAL30"  --filter-expression "QUAL < 30.0" '
        ' --filter-name "SOR3"  --filter-expression "SOR > 3.0" '
        ' --filter-name "FS60"  --filter-expression "FS > 60.0" '
        ' --filter-name "DP10"  --filter-expression "DP < 10.0" '
        ' --filter-name "MQRankSum-12.5"  --filter-expression "MQRankSum < -12.5" '
        ' --filter-name "ReadPosRankSum-8"  --filter-expression "ReadPosRankSum < -8.0" '





##################################







### Standardized Variant Filtering & Processing  ###

##### NOTE: These will be the same across ALL ALIGNER-VariantCaller Combos! #####


### Filter out all INDELs with length greater than 15 bp 
rule Filter_VarCaller_VCF_RemoveIndelsGreaterThan15bp:
    input:
        VarCaller_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.FilterTagged.vcf",
    output:
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bcftools view --types snps,indels -i 'abs(strlen(ALT)-strlen(REF))<=15' -f PASS,. {input.VarCaller_VCF} > {output.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only}"



# Perform masking of PLC regions

rule RemoveCoscollaRegions_From_VarCaller_VCF_RemoveIndelsGreaterThan15bp:
    input:
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        Coscolla_etal_Regions_BED = "References/Mtb_H37Rv_MaskingSchemes/201027_Mtb_H37rv_pLC_Regions_CoscollaExcludedGenes.bed"
    output:
        VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bedtools intersect -header -v -a {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -b {input.Coscolla_etal_Regions_BED} -wa > {output.VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} "


# Perform masking of Pmap-K50E4 < 1 (Pileup Mappability) regions

rule Remove_Pmap_K50E4_Below1_From_VarCaller_VCF_RemoveIndelsGreaterThan15bp:
    input:
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        Pmap_K50E4_Below1_Regions_BED = "References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed",
    output:
        VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_WiIll}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
    conda:
        "CondaEnvs/Illumina_AlnAndVC_Benchmarking_V3_Conda.yml"
    shell:
        "bedtools intersect -header -v -a {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -b {input.Pmap_K50E4_Below1_Regions_BED} -wa > {output.VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} "








##### Running Hap.py Benchmarking - For ALL Aligner-VariantCaller combinations ######

##### NOTE: These rules will be the same across ALL ALIGNER-VariantCaller combinations! #####

# For running Hap.py we need to identify the quality thresholds that will be used 


### A) For Pilon, we will use MQ
### B) For Bcftools-call, we will use QUAL
### C) For VarScan, we will use ____
### D) For HaplotypeCaller, we will use CNN_2D




rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved_MQforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ.summary.csv",
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "





rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_PLCRegionsRemoved_MQforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads:
        1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "






rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed_MQforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "



rule VarCall_VCF_Happy_All_WiStratificationBy_Pmap_and_SV_4sets_V2_MQforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.MQforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.MQ --threads {threads} "
        " --stratification {input.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} --roc-regions HighPmap_NoSV,LowPmap_NoSV,LowPmap_WiSV,HighPmap_WiSV -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "

















rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved_QUALforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ.summary.csv",
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc QUAL --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "






rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_PLCRegionsRemoved_QUALforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads:
        1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc QUAL --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "






rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed_QUALforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc QUAL --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "



rule VarCall_VCF_Happy_All_WiStratificationBy_Pmap_and_SV_4sets_V2_QUALforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.QUALforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc QUAL --threads {threads} "
        " --stratification {input.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} --roc-regions HighPmap_NoSV,LowPmap_NoSV,LowPmap_WiSV,HighPmap_WiSV -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "

















# CNN_2DforQQ



rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved_CNN_2DforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ.summary.csv",
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.CNN_2D --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "








rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_PLCRegionsRemoved_CNN_2DforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads:
        1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.CNN_2D --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "






rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed_CNN_2DforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.CNN_2D --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "



rule VarCall_VCF_Happy_All_WiStratificationBy_Pmap_and_SV_4sets_V2_CNN_2DforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.CNN_2DforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.CNN_2D --threads {threads} "
        " --stratification {input.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} --roc-regions HighPmap_NoSV,LowPmap_NoSV,LowPmap_WiSV,HighPmap_WiSV -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "











# QD (Quality by Depth metric from GATK-HaplotypeCaller)
### Relevant Documentation from GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360051304751-QualByDepth
### QD should be used instead of QUAL 

rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved_QDforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ.summary.csv",
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.QD --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "








rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_PLCRegionsRemoved_QDforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.CoscollaRegionsRemoved.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads:
        1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_PLCRegionsRemoved/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_CoscollaRegionsRemoved} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.QD --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "






rule VarCall_VCF_Happy_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed_QDforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.Pmap_K50E4_Below1_Removed.vcf",
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_Pmap_K50E4_Below1Removed/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_SNPsINDELs_1to15bp_Pmap_K50E4_Below1Removed} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.QD --threads {threads} -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "



rule VarCall_VCF_Happy_All_WiStratificationBy_Pmap_and_SV_4sets_V2_QDforQQ:
    input:
        H37rv_FA = refGenome_FA_PATH,
        VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/VarCall_{variantCaller_ID}/{sampleID_Wi_Ill_And_PB}.IllPE.{aligner_ID}.{variantCaller_ID}.PASS.SNPsINDELs.Lengths_1to15bp.vcf",   
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ.summary.csv"
    conda:
        "CondaEnvs/happy_3_12.yml"
    threads: 1
    params:
        output_Prefix = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/AlignerAndCaller_Benchmarking/{aligner_ID}/Happy_Benchmarking_{aligner_ID}_{variantCaller_ID}_ShortVariants_NoRegionsRemoved_StratifiedBy_PmapAndSVs_4sets_SV_LengthFiltered_50bp/Hap.py.{sampleID_Wi_Ill_And_PB}.QDforQQ"
    shell:
        "export HGREF=/n/data1/hms/dbmi/farhat/mm774/References/Happy_Ref_hg19/hg19.fa \n"
        "export PATH=/n/data1/hms/dbmi/farhat/mm774/ProcessedData_Etc/Happy_VC_Eval_TestDir/rtg-tools-3.11/:$PATH \n"
        " "
        "hap.py {input.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} {input.VarCaller_VCF_Pass_SNPsAndINDELs_Lengths_1to15bp_Only} -r {input.H37rv_FA} " 
        "-o {params.output_Prefix} "
        "--engine=vcfeval --preprocess-truth  --pass-only --roc INFO.QD --threads {threads} "
        " --stratification {input.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} --roc-regions HighPmap_NoSV,LowPmap_NoSV,LowPmap_WiSV,HighPmap_WiSV -f {input.MM2_AtoRef_PP_NOT_AMB_Regions_BED} "













