
##################################################################################################################
############### Coverage Analysis Pre-processing - PBassembly_Flye_I3_PP_CoverageAnalysis ########################
##################################################################################################################


# 1 #
############ MM2: Align PacBio Subreads To PacBio Flye I3 PP Assembly ################################


rule align_PacBio_Subreads_To_Flye_I3_PP_Assembly_With_Minimap2:
    input:
        pb_subreads_fq = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/pacbio_reads/{sampleID_Wi_Ill_And_PB}.merged.subreads.fastq.gz",
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
    output:
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_SAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.sam",
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam",
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam.bai",
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM_Depth_TXT = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam.depth.txt"

    threads: 10
    conda:
        "CondaEnvs/PacBio_Software_py27_Conda.yml"
    shell:
        "minimap2 -t {threads} -ax map-pb {input.I3_Assembly_PilonPolished_FA} {input.pb_subreads_fq} > {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_SAM} \n"
        "samtools view -bS {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_SAM} | samtools sort - > {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM} \n"
        "samtools index {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM} \n"
        "samtools depth -a {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM} > {output.pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM_Depth_TXT}"



rule pilon_VarCalling_PacBio_Subreads_AlignTo_Flye_I3_PP_Assembly_DefaultParam:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam",
        pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/PacBio_Subreads_AlignedTo_Flye_I3_PP_Minimap2/{sampleID_Wi_Ill_And_PB}.pb.subreads.AlnTo.Flye_I3_PP.minimap2.bam.bai",
    output:
        pilon_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_PacBio_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/{sampleID_Wi_Ill_And_PB}.PacBio.AlnToFlye_I3_PP.vcf",
        pilon_ChangesFile = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_PacBio_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/{sampleID_Wi_Ill_And_PB}.PacBio.AlnToFlye_I3_PP.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_PacBio_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx35g --genome {input.I3_Assembly_PilonPolished_FA} --pacbio  {input.pb_subreads_MM2_To_Flye_I3_PP_Assembly_BAM} --output {wildcards.sampleID_Wi_Ill_And_PB}.PacBio.AlnToFlye_I3_PP"
        " --outdir {params.Pilon_OutputDir_PATH} --changes --tracks --fix snps,indels --vcf"








# 2 #
############ BWA-MEM: Align Illumina WGS reads To PacBio Flye I3 PP Assembly ################################


rule bwa_idx_ref_Flye_I3_PP_Assembly:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
    output:
        I3_Assembly_PilonPolished_FA_BWT = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta.bwt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "bwa index {input}"


rule bwa_map_IllPE_AlignTo_Flye_I3_PP_Assembly:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        I3_Assembly_PilonPolished_FA_BWT = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta.bwt",
        fq1_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.sam"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_Wi_Ill_And_PB}\tSM:{sampleID_Wi_Ill_And_PB}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.I3_Assembly_PilonPolished_FA} {input.fq1_trimmed} {input.fq2_trimmed} > {output}"


rule samtools_ViewSortAndGetDepth_IllPE_AlignTo_Flye_I3_PP:
    input:
        Illumina_AlignTo_Flye_I3_PP_SAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.sam",
    output:
        Illumina_AlignTo_Flye_I3_PP_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.bam",
        Illumina_AlignTo_Flye_I3_PP_BAM_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.bam.bai",
        Illumina_AlignTo_Flye_I3_PP_BAM_Depth_TXT = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.bam.depth.txt",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input.Illumina_AlignTo_Flye_I3_PP_SAM} | samtools sort - > {output.Illumina_AlignTo_Flye_I3_PP_BAM} \n"
        "samtools index {output.Illumina_AlignTo_Flye_I3_PP_BAM} \n"
        "samtools depth -a {output.Illumina_AlignTo_Flye_I3_PP_BAM} > {output.Illumina_AlignTo_Flye_I3_PP_BAM_Depth_TXT}"


#### PICARD (remove duplicates) AND Get Depth with Samtools #####

rule picard_RemoveDup_AND_samtools_ViewSortAndGetDepth_IllPE_AlignTo_Flye_I3_PP:
    input:
        Illumina_AlignTo_Flye_I3_PP_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.bam",
    output:
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_METRICS = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.metrics",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.bai",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_Depth_TXT = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.depth.txt",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_Depth_Average_TXT = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.depth.averaged.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx4g MarkDuplicates I={input.Illumina_AlignTo_Flye_I3_PP_BAM} O={output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_METRICS} ASSUME_SORT_ORDER=coordinate \n"
        " "
        "samtools index {output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM} \n"
        "samtools depth -a {output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM} > {output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_Depth_TXT} \n"
        "samtools depth -a {output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM} | awk '{{sum+=$3}} END {{ print \"Average = \",sum/NR}}' > {output.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_Depth_Average_TXT}"





rule pilon_VarCalling_IllPE_AlignTo_AlignTo_Flye_I3_PP_DefaultVariantCalling:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam",
        Illumina_AlignTo_Flye_I3_PP_Duprem_BAM_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/IlluminaPE_AlignedTo_Flye_I3_PP_bwamem/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.Flye_I3_PP.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_IlluminaPE_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/{sampleID_Wi_Ill_And_PB}.IllPE.FlyeI3_Assembly.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_IlluminaPE_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/{sampleID_Wi_Ill_And_PB}.IllPE.FlyeI3_Assembly.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_Wi_Ill_And_PB}/PBassembly_Flye_I3_PP_CoverageAnalysis/Pilon_IlluminaPE_AlignedTo_FlyeI3_Assembly_DefaultParameters_VariantCalling/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx14g --genome {input.I3_Assembly_PilonPolished_FA} --bam {input.Illumina_AlignTo_Flye_I3_PP_Duprem_BAM} --output {wildcards.sampleID_Wi_Ill_And_PB}.IllPE.FlyeI3_Assembly"
        " --outdir {params.Pilon_OutputDir_PATH} --changes --tracks  --variant"





