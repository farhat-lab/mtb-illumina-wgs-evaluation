

###### Illumina Analysis  #####


rule unzip_FQ_GZs:
    input:
        fq1_gz = lambda wildcards: SampleIDTo_Illumina_PE_FQ1_Dict[wildcards.sampleID_WiIll],
        fq2_gz = lambda wildcards: SampleIDTo_Illumina_PE_FQ2_Dict[wildcards.sampleID_WiIll],
    output:
        fq1_unzipped = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq",
        fq2_unzipped = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq",
    threads: 1
    shell:
        "zcat {input.fq1_gz} > {output.fq1_unzipped} \n"
        "zcat {input.fq2_gz} > {output.fq2_unzipped} \n"


# Define read number for forward and reverse reads
readNum = ['1', '2']

rule fastqc:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_{readNum}.fastq",
    output:
        html = output_Dir + "{sampleID_WiIll}/IlluminaWGS/fastqc/{sampleID_WiIll}_fq{readNum}_fastqc.html",
        zip = output_Dir + "{sampleID_WiIll}/IlluminaWGS/fastqc/{sampleID_WiIll}_fq{readNum}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        output_Dir+"logs/fastqc/{sampleID_WiIll}_fq{readNum}.log"
    wrapper:
        "0.38.0/bio/fastqc"

rule multiqc_AllDirs:
    input:
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/fastqc/{sampleID_WiIll}_fq{readNum}_fastqc.zip", sampleID_WiIll=input_SampleIDs_WithIllumina, readNum=['1', '2']),
    output:
        output_Dir + "/multiqc_Reports/multiqc_AllDirs.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        output_Dir + "/logs/multiqc/multiqc_AllDirs.log"
    wrapper:
        "0.38.0/bio/multiqc"






# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_PE_Trimming:
    input:
        r1 = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq",
        r2 = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq",
    output:
        r1 = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        r2 = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.unpaired.fastq",
        r2_unpaired = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.unpaired.fastq",
    log:
        output_Dir + "/logs/trimmomatic/{sampleID_WiIll}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'References/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        # optional parameters
        # extra=" "
    threads: 8
    wrapper:
        "0.38.0/bio/trimmomatic/pe"



# Define read number for forward and reverse reads
readNum = ['1', '2']

rule fastqc_OfTrimmomaticReads:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_{readNum}_trimmed.fastq",
    output:
        html = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/fastqc_Trimmomatic/{sampleID_WiIll}_fq{readNum}_fastqc.html",
        zip = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/fastqc_Trimmomatic/{sampleID_WiIll}_fq{readNum}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        output_Dir +"/logs/fastqc_Trimmomatic/{sampleID_WiIll}_fq{readNum}.log"
    wrapper:
        "0.38.0/bio/fastqc"


rule multiqc_OfTrimmomaticReads:
    input:
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/fastqc_Trimmomatic/{sampleID_WiIll}_fq{readNum}_fastqc.zip", sampleID_WiIll=input_SampleIDs_WithIllumina, readNum=['1', '2']),
        expand(output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.stats.txt", sampleID_WiIll=input_SampleIDs_WithIllumina),

    output:
        output_Dir + "multiqc_Reports/multiqc_Trimmomatic.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        output_Dir+"/logs/multiqc/multiqc_Trimmomatic.log"
    wrapper:
        "0.38.0/bio/multiqc"


rule bwa_map_IllPE_AlignTo_H37rv:
    input:
        fa = refGenome_FA_PATH,
        fq1_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
        r1_unpaired = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.unpaired.fastq",
        r2_unpaired = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.unpaired.fastq",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.sam"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_WiIll}\tSM:{sampleID_WiIll}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.fa} {input.fq1_trimmed} {input.fq2_trimmed} > {output}"




rule samtools_ViewAndSort_IllPE_AlignTo_H37rv:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.sam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam"
    log:
        output_Dir + "/logs/IlluminaWGS/samtools_sort_IllPE_AlignTo_H37rv/{sampleID_WiIll}.samtools.sort.log"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} -m 4G | samtools sort -m 4G - > {output}"



rule samtools_index_IllPE_AlignTo_H37rv:
    input:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam"
    output:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam.bai"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"


rule samtools_Depth_AverageAll_IllPE_AlignTo_H37rv:
    input:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam"
    output:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam.depth.averaged.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools depth -a {input} -m 8G | awk '{{sum+=$3}} END {{ print \"Average = \",sum/NR}}' > {output}"


# "awk '$4 == \"+\" {{print $1}}' {input.assembly_info_txt} > {output.assembly_circularcontigs_txt} \n"



#####################################
#### PICARD (remove duplicates) #####
#####################################

rule picard_RemoveDup_IllPE_AlignTo_H37rv:
    input:
        IllPE_BwaMEM_BAM = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.bam",
    output:
        IllPE_BwaMEM_Duprem_BAM = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
        IllPE_BwaMEM_Duprem_METRICS = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.metrics",

    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx6g MarkDuplicates I={input.IllPE_BwaMEM_BAM} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"



rule samtools_index_IllPE_Duprem_AlignTo_H37rv:
    input:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
    output:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.bai",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"



rule samtools_Depth_AverageAll_IllPE_AlignTo_H37rv_Duprem:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam"
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.depth.averaged.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools depth -a {input} | awk '{{sum+=$3}} END {{ print \"Average = \",sum/NR}}' > {output}"


rule samtools_Depth_IllPE_AlignTo_H37rv_Duprem:
    input:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
    output:
        output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.depth.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools depth -a {input} > {output}"

rule samtools_Stats_IllPE_AlignTo_H37rv_Duprem:
    input:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
    output:
        output_Dir + "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.stats.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools stats {input} > {output}"







rule pilon_VarCalling_IllPE_AlignTo_H37rv_Default_Variant:
    input:
        Ref_fa = refGenome_FA_PATH,
        IllPE_BwaMEM_Duprem_BAM = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
        IllPE_BwaMEM_Duprem_BAM_BAI = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiIll}.IllPE.H37rv.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiIll}.IllPE.H37rv.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_DefaultParameters_VariantCalling/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx12g --genome {input.Ref_fa} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_WiIll}.IllPE.H37rv"
        " --outdir {params.Pilon_OutputDir_PATH} --changes --tracks  --variant"


# Using Pilon with "--fix all,breaks --minmq 1 --mindepth 5" 

rule pilon_VarCalling_IllPE_AlignTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks:
    input:
        Ref_fa = refGenome_FA_PATH,
        IllPE_BwaMEM_Duprem_BAM = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
        IllPE_BwaMEM_Duprem_BAM_BAI = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx12g --genome {input.Ref_fa} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_WiIll}.IllPE.H37rv"
        " --outdir {params.Pilon_OutputDir_PATH} --fix all,breaks --vcf --changes --tracks --minmq 1 --mindepth 5"



rule Filter_Pilon_VarCalling_VCF_ToH37rv:
    input:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf",
    output:
        pilon_VCF_PassOnly = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.AllVariants.PassOnly.vcf",
        pilon_VCF_SNPsOnly_PassOnly = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.PassOnly.SNPsOnly.vcf",
        pilon_VCF_SNPsOnly = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.SNPsOnly.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view --types snps,indels,mnps,other -f PASS {input.pilon_VCF} > {output.pilon_VCF_PassOnly} \n"
        "bcftools view --types snps -f PASS {input.pilon_VCF} > {output.pilon_VCF_SNPsOnly_PassOnly} \n"
        "bcftools view --types snps {input.pilon_VCF} > {output.pilon_VCF_SNPsOnly} \n"



rule bcftools_mpileup_VarCalling_IllPE_AlignTo_H37rv:
    input:
        H37rv_FA = refGenome_FA_PATH,
        IllPE_BwaMEM_Duprem_BAM = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam",
        IllPE_BwaMEM_Duprem_BAM_BAI = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.IllPE.H37rv.duprem.bam.bai",
    output:
        mpileup_VCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/bcftools_mpileup_call_IlluminaPE_AlignedTo_H37rv/{sampleID_WiIll}.bcftools.call.IllPE.H37rv.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "samtools mpileup -ugf {input.H37rv_FA} {input.IllPE_BwaMEM_Duprem_BAM}"
        " | bcftools call -vm -O v -o {output.mpileup_VCF} "






rule calculate_F2_Score_IlluminaPE:
    input:
        pilon_VCF = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.vcf",
        lineage_def_ref_pos_TXT = "Scripts/positions_to_select.txt",
        Coll2014_LinSpeSNPs_final_CSV = "Scripts/Coll2014_LinSpeSNPs_final.csv",        
    output:
        F2_Score_TXT = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/{sampleID_WiIll}_F2.txt",
        #TEMP_BCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf",
        #TEMP_BCF_CSI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf",
    params:
        F2_OutputDir_PATH = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation",
        TEMP_BCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf",
        TEMP_BCF_CSI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf.csi",
    shell:
        "python Scripts/lineage_defining_SNP_depths_collection_and_F2_calculation.py {input.pilon_VCF} {params.F2_OutputDir_PATH} {input.lineage_def_ref_pos_TXT} {wildcards.sampleID_WiIll} {input.Coll2014_LinSpeSNPs_final_CSV} \n"
        "rm {params.TEMP_BCF} {params.TEMP_BCF_CSI} "






########################################################
### LINEAGE CALLING script using VCF as input ###
########################################################

rule LineageCall_Coll_IlluminaWGS_AlignTo_H37rv:
    input:
        pilon_VCF_SNPsOnly_PassOnly = output_Dir+ "/{sampleID_WiIll}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_WiIll}.IllPE.H37rv.PassOnly.SNPsOnly.vcf",
        Coll_SNPs_PATH = "References/LineageCalling_SNP_Schemes/coll.tsv"
    output:
        output_Dir + "/{sampleID_WiIll}/LineageCalling/LineageCall_IlluminaWGS_AlignTo_H37rv/{sampleID_WiIll}.IlluminaWGS.Pilon.lineage_call.Coll.txt"
    shell:
        "Scripts/fast-lineage-caller-vcf.py {input.pilon_VCF_SNPsOnly_PassOnly} {input.Coll_SNPs_PATH} > {output} "




