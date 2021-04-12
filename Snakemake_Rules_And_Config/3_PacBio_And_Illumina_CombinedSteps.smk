




##### Combined Illumina + PacBio Analysis Steps #####

rule bwa_idx_ref_100kbContigs_FlyeAssembly_I3:
    input:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta.bwt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "bwa index {input}"


rule bwa_map_IllPE_AlignTo_I3_Assembly:
    input:
        PacBio_Flye_Assembly_I3_FA = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
        PacBio_Flye_Assembly_I3_FA_BWT_IDX = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta.bwt",
        fq1_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.sam"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_Wi_Ill_And_PB}\tSM:{sampleID_Wi_Ill_And_PB}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.PacBio_Flye_Assembly_I3_FA} "
        "{input.fq1_trimmed} {input.fq2_trimmed} > {output}"


rule samtools_ViewAndSort_IllPE_AlignTo_I3_Assembly:
    input:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.sam"
    output:
        bam = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam",
        bai = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam.bai",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} | samtools sort - > {output.bam} \n"
        "samtools index {output.bam}"



#####################################
#### PICARD (remove duplicates) #####
#####################################

rule picard_RemoveDup_IllPE_AlignTo_I3_Assembly:
    input:
        IllPE_BwaMEM_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam",
    output:
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.bai",
        IllPE_BwaMEM_Duprem_METRICS = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.metrics",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx5g MarkDuplicates I={input.IllPE_BwaMEM_BAM} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"
        " \n"
        "samtools index {output.IllPE_BwaMEM_Duprem_BAM}"


rule pilon_IllPE_Polishing_I3_Assembly:
    input:
        PacBio_Flye_Assembly_I3_FA = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.vcf",
        pilon_I3_Polished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        pilon_I3_PP_ChangesFile = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.changes"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx10g --fix snps,indels --genome {input.PacBio_Flye_Assembly_I3_FA} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished"
        " --outdir {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/ --variant --changes --tracks \n"
        "samtools faidx {output.pilon_I3_Polished_FA}"



rule samtools_faidx_PilonPolished_I3Assembly:
    input:
        output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
    output:
        output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta.fai"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    threads: 1
    shell:
        "samtools faidx {input}"


rule FastANI_I3_PP_Assembly: 
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        H37rv_Ref_fa = refGenome_FA_PATH,     
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/FastANI_OutputDirs/Flye_I3_PilonPolished_FastANI/{sampleID_Wi_Ill_And_PB}_I3_PP_AssemblyToH37rv.FastANI.txt"
    conda:
        "CondaEnvs/fastani_1_3_0_Conda.yml"
    shell:
        "fastANI -q {input.I3_Assembly_PilonPolished_FA} -r {input.H37rv_Ref_fa} -o {output}"


rule Prokka_Anno_I3_PilonPolished:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
    output:
        I3_Prokka_Anno_FNA = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Prokka/Flye_I3_PB_PilonPolished_Prokka_Anno/{sampleID_Wi_Ill_And_PB}.fna"
    conda:
        "CondaEnvs/Prokka_1130_Conda.yml" #"CondaEnvs/Prokka_1145_Conda.yml"
    threads: 8
    shell: 
        "prokka --force --cpus {threads} --genus Mycobacterium --species tuberculosis --strain {wildcards.sampleID_Wi_Ill_And_PB}"
        " --addgenes --mincontiglen 100000 --centre XXX --proteins {H37rv_GBK_PATH} "
        " --outdir {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/Prokka/Flye_I3_PB_PilonPolished_Prokka_Anno/ --prefix {wildcards.sampleID_Wi_Ill_And_PB} {input.I3_Assembly_PilonPolished_FA}"













#########################################################
#########################################################
## Analysis Versus H37rv (call variants against H37Rv) ##
#########################################################
#########################################################


##############################################################################
############ Minimap2: Assembly To H37rv Alignment & Variant Calling #########
##############################################################################

rule Minimap2_Flye_I3_PP_AlignTo_H37rv:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_Flye_I3_PP_To_H37rv_SAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.sam",
        MM2_Flye_I3_PP_To_H37rv_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam",
        MM2_Flye_I3_PP_To_H37rv_bai = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.bai",
        MM2_Flye_I3_PP_To_H37rv_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.vcf",
        MM2_Flye_I3_PP_To_H37rv_PAF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paf",
    conda:
        "CondaEnvs/PacBio_Software_py27_Conda.yml"
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
    shell: 
        "minimap2 -ax asm10 --cs {input.H37rv_FA} {input.I3_Assembly_PilonPolished_FA} | awk '$1 ~ /^@/ || ($5 == 60)' > {output.MM2_Flye_I3_PP_To_H37rv_SAM} \n"
        "samtools view -bS {output.MM2_Flye_I3_PP_To_H37rv_SAM} | samtools sort - > {output.MM2_Flye_I3_PP_To_H37rv_BAM} \n"
        "samtools index {output.MM2_Flye_I3_PP_To_H37rv_BAM} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.I3_Assembly_PilonPolished_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_Wi_Ill_And_PB} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_Flye_I3_PP_To_H37rv_VCF} \n"                         
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.I3_Assembly_PilonPolished_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_Wi_Ill_And_PB} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} - > {output.MM2_Flye_I3_PP_To_H37rv_PAF} \n"




rule filter_MM2_paftools_VCF_GC3_PP_AlignTo_H37rv_RemoveIndelsGreaterThan15bp:
    input:
        MM2_Flye_I3_To_H37rv_paftools_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.vcf",
    output:
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_Flye_I3_To_H37rv_paftools_VCF_SNPs_Only = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.SNPsOnly.vcf",    
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell:
        "bcftools view --types snps,indels -i 'abs(strlen(ALT)-strlen(REF))<=15' {input.MM2_Flye_I3_To_H37rv_paftools_VCF} > {output.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} \n"
        "bcftools view --types snps {input.MM2_Flye_I3_To_H37rv_paftools_VCF} > {output.MM2_Flye_I3_To_H37rv_paftools_VCF_SNPs_Only} \n"




rule bcftools_mpileup_MM2BAM_GC3_PP_AlignTo_H37rv:
    input:
        MM2_Flye_I3_PP_To_H37rv_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam",
        H37rv_FA = refGenome_FA_PATH,     
    output:
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_out = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf",
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.KeepAllPositions.vcf"
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell: # Run bcftools mpileup to summarize coverage and basepair outputs from the Minimap2 alignment (BAM)
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_Flye_I3_PP_To_H37rv_BAM} > {output.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_out} \n"
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_Flye_I3_PP_To_H37rv_BAM} | bcftools call -c -o {output.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF}"




tab_str="\t"

rule infer_DuplicatedRegions_from_mpileup_MM2BAM_GC3_PP_AlignTo_H37rv:
    input:
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_out = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf",
        H37rv_bedtools_genome="References/bedtools_ref/GCF_000195955.2_ASM19595v2_genomic.fasta.bedtools.genome",

    output:
        MM2_AtoRef_PP_DupRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf.DupRegions.bed",
        MM2_AtoRef_PP_NOTDupRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf.NOTDupRegions.bed",
        DupRegions_stratification_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}_DUP_Regions_stratification.tsv"
    threads: 1
    shell: 
        "bcftools view -i 'DP>=2' {input.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_out} | vcf2bed | bedtools merge > {output.MM2_AtoRef_PP_DupRegions_BED} \n"
        "bedtools complement -i {output.MM2_AtoRef_PP_DupRegions_BED} -g {input.H37rv_bedtools_genome} > {output.MM2_AtoRef_PP_NOTDupRegions_BED} \n"
        "echo -e 'Duplicated{tab_str}{wildcards.sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf.DupRegions.bed' > {output.DupRegions_stratification_TSV}   \n"
        "echo -e 'NotDuplicated{tab_str}{wildcards.sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf.NOTDupRegions.bed' >> {output.DupRegions_stratification_TSV}   \n"










rule Minimap2_Flye_I3_PBonly_AlignTo_H37rv:
    input:
        PacBio_Flye_Assembly_I3_FA = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta",
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_Flye_I3_PBonly_To_H37rv_SAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.sam",
        MM2_Flye_I3_PBonly_To_H37rv_BAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.bam",
        MM2_Flye_I3_PBonly_To_H37rv_bai = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.bam.bai",
        MM2_Flye_I3_PBonly_To_H37rv_VCF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.vcf",
        MM2_Flye_I3_PBonly_To_H37rv_PAF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.paf",
    conda:
        "CondaEnvs/PacBio_Software_py27_Conda.yml"
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
    shell: 
        "minimap2 -ax asm10 --cs {input.H37rv_FA} {input.PacBio_Flye_Assembly_I3_FA} > {output.MM2_Flye_I3_PBonly_To_H37rv_SAM} \n"
        "samtools view -bS {output.MM2_Flye_I3_PBonly_To_H37rv_SAM} | samtools sort - > {output.MM2_Flye_I3_PBonly_To_H37rv_BAM} \n"
        "samtools index {output.MM2_Flye_I3_PBonly_To_H37rv_BAM} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.PacBio_Flye_Assembly_I3_FA} | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiPacBio} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_Flye_I3_PBonly_To_H37rv_VCF} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.PacBio_Flye_Assembly_I3_FA} | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiPacBio} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} - > {output.MM2_Flye_I3_PBonly_To_H37rv_PAF} \n"




rule Filter_SNPsOnly_MM2_VCF_GC3_PBonly_ToH37rv:
    input:
        MM2_Flye_I3_PBonly_To_H37rv_VCF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.vcf",
    output:
        MM2_Flye_I3_PBonly_To_H37rv_VCF_SNPsOnly = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.SNPsOnly.vcf",
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    shell:
        "bcftools view --types snps {input.MM2_Flye_I3_PBonly_To_H37rv_VCF} > {output.MM2_Flye_I3_PBonly_To_H37rv_VCF_SNPsOnly} "




rule filter_MM2_paftools_VCF_I3_PBonly_AlignTo_H37rv_RemoveIndelsGreaterThan15bp:
    input:
        MM2_Flye_I3_PBonly_To_H37rv_VCF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.vcf",
    output:
        MM2_Flye_I3_PBonly_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.paftools.Lengths_1to15bp.vcf",    
        MM2_Flye_I3_PBonly_To_H37rv_paftools_VCF_SNPs_Only = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.paftools.SNPsOnly.vcf",    
    conda:
        "CondaEnvs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell:
        "bcftools view --types snps,indels -i 'abs(strlen(ALT)-strlen(REF))<=15' {input.MM2_Flye_I3_PBonly_To_H37rv_VCF} > {output.MM2_Flye_I3_PBonly_To_H37rv_paftools_VCF_SNPsAndINDELs_Lengths_1to15bp_Only} \n"
        "bcftools view --types snps {input.MM2_Flye_I3_PBonly_To_H37rv_VCF} > {output.MM2_Flye_I3_PBonly_To_H37rv_paftools_VCF_SNPs_Only} \n"

















bioawk_PadBy100bp_string="'{$start=$start - 100 ; $end=$end + 100 ;  print $0 }'"

rule NucDiff_Analysis_Flye_I3_PP_vs_H37rv_V2:
    input:
        I3_Assembly_PilonPolished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        H37rv_FA = refGenome_FA_PATH,
        H37rv_bedtools_genome="References/bedtools_ref/GCF_000195955.2_ASM19595v2_genomic.fasta.bedtools.genome",
        #Q3Assembly_PilonPolished_FA = output_Dir + "{sampleID_Wi_Ill_And_PB}/GC3_IlluminaPolishing/pilon_IllPE_Polishing_Q3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.IllPE.Q3Assembly.fasta",
    output:
        NucDiff_GC3_PP_To_H37rv_SNPs_GFF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_snps.gff",
        NucDiff_GC3_PP_To_H37rv_SNPs_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_snps.vcf",
        NucDiff_GC3_PP_To_H37rv_Struct_GFF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.gff",
        NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_GFF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.gff",
        NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_DeletionsOnly_GFF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.DeletionsOnly.gff",
        NucDiff_GC3_PP_To_H37rv_Struct_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.bed", 
        NucDiff_FilteredSVs_With50bpLengthFilter_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.bed", 
        NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.NoInversions.bed",
        WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.NoInversions.PaddedBy100bp.bed",
        NoSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_COMPLEMENT_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.NoInversions.PaddedBy100bp.complement.bed",
    conda:  
        "CondaEnvs/nucdiff_2_0_3_Conda.yml"
    threads:
        1
    shell:
        "nucdiff {input.H37rv_FA} {input.I3_Assembly_PilonPolished_FA} {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{wildcards.sampleID_Wi_Ill_And_PB} NucDiff_{wildcards.sampleID_Wi_Ill_And_PB} --vcf yes \n"
        " " # Filter out reshuffling SVs from NucDiff's output
        ' grep -v "reshuffling" {output.NucDiff_GC3_PP_To_H37rv_Struct_GFF} > {output.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_GFF} \n'


        " " # Use custom script to convert gff to bed format
        "python Scripts/Parse_NucDiff_gff_to_bed.py "
        "--nucdiff_sv_gff {output.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_GFF} " # arg1
        " --sample_id {wildcards.sampleID_Wi_Ill_And_PB} \n"

        " " # use bioawk to pad the start and end of each BED entry by 100 bp
        "bioawk -c bed {bioawk_PadBy100bp_string} {output.NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_BED} | sed -e 's/-100/0/g' > {output.WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED} \n"
        

        " " # Get complement of NucDiff SV regions
        "bedtools complement -i {output.WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED} -g {input.H37rv_bedtools_genome} > {output.NoSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_COMPLEMENT_BED} \n"                                                                                                                             
        
        " "
        " " # Filter NucDiff SVs for DELETIONS ONLY
        "grep -E 'deletion|collapsed_repeat|collapsed_tandem_repeat' {output.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_GFF} > {output.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_DeletionsOnly_GFF}"





        








tab_str="\t"

rule process_Pmap_and_NucDiffSV_Regions_ForStratification_V2_With50bpLengthFilter_SVs:
    input:
        H37rv_bedtools_genome="References/bedtools_ref/GCF_000195955.2_ASM19595v2_genomic.fasta.bedtools.genome",
        Pmap_K50E4_Below1_Regions_BED = "References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed",
        WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.NoInversions.PaddedBy100bp.bed",
        NoSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_COMPLEMENT_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.Parsed.With50bpLengthFilter.NoInversions.PaddedBy100bp.complement.bed",
    output:
        Low_PileupMappability_Regions = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/Low_PileupMappabilityRegions.PMap_K50E4_BELOW_1.bed",
        High_PileupMappability_Regions = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/High_PileupMappabilityRegions.PMap_K50E4_ABOVE_1.bed",
        i_HighMap_NoSV_IntersectionOfRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/StratRegions.V2.HighMap_NoSV_Regions.intersected.{sampleID_Wi_Ill_And_PB}.bed",
        i_LowMap_NoSV_IntersectionOfRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/StratRegions.V2.LowMap_NoSV_Regions.intersected.{sampleID_Wi_Ill_And_PB}.bed",
        i_LowMap_WiSV_IntersectionOfRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/StratRegions.V2.LowMap_WiSV_Regions.intersected.{sampleID_Wi_Ill_And_PB}.bed",
        i_HighMap_WiSV_IntersectionOfRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/StratRegions.V2.HighMap_WiSV_Regions.intersected.{sampleID_Wi_Ill_And_PB}.bed",
        i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/Hap.py_VariantCalling_EvalDir/Happy_StratificationFiles/stratificationRegions.V2.{sampleID_Wi_Ill_And_PB}.tsv",
    conda:
        "CondaEnvs/nucdiff_2_0_3_Conda.yml"
    threads: 1
    shell:
        " " # Copy and complement the Low Pileup Mappability regions
        "cp {input.Pmap_K50E4_Below1_Regions_BED} {output.Low_PileupMappability_Regions} \n"
        " " # Get complement of low Mappability regions (Thus generating the HIGH mappability regions)
        "bedtools complement -i {output.Low_PileupMappability_Regions} -g {input.H37rv_bedtools_genome} > {output.High_PileupMappability_Regions} \n"
        "bedtools intersect -a {output.High_PileupMappability_Regions} -b {input.NoSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_COMPLEMENT_BED} > {output.i_HighMap_NoSV_IntersectionOfRegions_BED} \n"
        "bedtools intersect -a {output.Low_PileupMappability_Regions} -b {input.NoSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_COMPLEMENT_BED} > {output.i_LowMap_NoSV_IntersectionOfRegions_BED} \n"
        "bedtools intersect -a {output.Low_PileupMappability_Regions} -b {input.WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED} > {output.i_LowMap_WiSV_IntersectionOfRegions_BED} \n"
        "bedtools intersect -a {output.High_PileupMappability_Regions} -b {input.WiSV_Regions_NucDiff_FilteredSVs_With50bpLengthFilter_NoInversions_PaddedBy100bp_BED} > {output.i_HighMap_WiSV_IntersectionOfRegions_BED} \n"
        "\n"
        "echo -e 'HighPmap_NoSV{tab_str}StratRegions.HighMap_NoSV_Regions.intersected.{wildcards.sampleID_Wi_Ill_And_PB}.bed' > {output.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} \n"
        "echo -e 'LowPmap_NoSV{tab_str}StratRegions.LowMap_NoSV_Regions.intersected.{wildcards.sampleID_Wi_Ill_And_PB}.bed' >> {output.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} \n"
        "echo -e 'LowPmap_WiSV{tab_str}StratRegions.LowMap_WiSV_Regions.intersected.{wildcards.sampleID_Wi_Ill_And_PB}.bed' >> {output.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} \n"
        "echo -e 'HighPmap_WiSV{tab_str}StratRegions.HighMap_WiSV_Regions.intersected.{wildcards.sampleID_Wi_Ill_And_PB}.bed' >> {output.i_stratRegions_V2_By_Pmap_and_SVs_ForHappy_TSV} \n"










 
rule calcEBR_EmpiricalBasePairRecall_V7:
    input:
        MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
        NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_DeletionsOnly_GFF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_Wi_Ill_And_PB}/results/NucDiff_{sampleID_Wi_Ill_And_PB}_ref_struct.Filtered.SVs.DeletionsOnly.gff",
        pilon_VCF = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/Pilon_IlluminaPE_AlignedTo_H37rv_minMQ_1_minDP_5_Fix_All_Breaks/{sampleID_Wi_Ill_And_PB}.IllPE.H37rv.vcf",
    output:
        EBR_Indiv_TSV = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.tsv",
        EBR_Indiv_OutcomeBreakdown_Dict_JSON = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.OutcomeBreakdown.json",
        EBR_Indiv_BED_Regions = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.bed",
        EBR_Indiv_BED_AmbOnlyRegions = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.Ambiguous.Regions.bed",  
    shell:
        "python Scripts/calculate_EBR_V7.py "
        "{input.MM2_AtoRef_Flye_I3_PP_BAM_mpileup_call_KeepAllPositions_VCF} " # arg1
        "{input.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_DeletionsOnly_GFF} "  # arg2
        "{input.pilon_VCF} " # arg3
        "{output.EBR_Indiv_TSV} {output.EBR_Indiv_OutcomeBreakdown_Dict_JSON} " # arg4 and arg5
        " {output.EBR_Indiv_BED_Regions} {output.EBR_Indiv_BED_AmbOnlyRegions}  " # arg6 and arg7





tab_str="\t"

rule infer_NonAmbiguous_from_Amb_EBR_Regions:
    input:
        EBR_Indiv_BED_AmbOnlyRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.Ambiguous.Regions.bed",
        H37rv_bedtools_genome="References/bedtools_ref/GCF_000195955.2_ASM19595v2_genomic.fasta.bedtools.genome",
    output:
        MM2_AtoRef_PP_NOT_AMB_Regions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.NOT.Ambiguous.Regions.bed",
    threads: 1
    shell: 
        "bedtools complement -i {input.EBR_Indiv_BED_AmbOnlyRegions_BED} -g {input.H37rv_bedtools_genome} > {output.MM2_AtoRef_PP_NOT_AMB_Regions_BED} \n"
 



rule BED_substract_DupRegions_From_AmbRegions:
    input:
        EBR_Indiv_BED_AmbOnlyRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.Ambiguous.Regions.bed",
        MM2_AtoRef_PP_DupRegions_BED = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.bam.mpileup.txt.vcf.DupRegions.bed",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/EmpiricalBasePairRecall_Analysis_V7_PacBio_Vs_IlluminaPilon/EBR.V7.IndivIsolate.{sampleID_Wi_Ill_And_PB}.Ambiguous.NoDupRegions.HighSeqDivergenceOnly.Regions.bed",
    threads: 1
    shell: 
       "bedtools subtract -a {input.EBR_Indiv_BED_AmbOnlyRegions_BED} -b {input.MM2_AtoRef_PP_DupRegions_BED} | cut -f 1-3 > {output}"






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



# Perform filtering of Pilon variant calls (1 to 15 bp Variants)

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







rule Happy_AmbRemoved_VCeval_T_PB_G3PP_MM2_Paftools_Vs_Q_Ill_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCs_SNPsINDELs_Lengths_1to15bp_NoRegionsRemoved:
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




rule Happy_AmbRemoved_VCeval_T_PB_G3PP_MM2_Paftools_Vs_Q_Ill_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCs_SNPsINDELs_Lengths_1to15bp_CoscollaRegionsRemoved:
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





rule Happy_AmbRemoved_VCeval_T_PB_G3PP_MM2_Paftools_Vs_Q_Ill_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCs_SNPsINDELs_Lengths_1to15bp_Pmap_K50E4_Below1_Removed:
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





rule Happy_AmbRemoved_VCeval_T_PB_G3PP_MM2_paftools_Vs_Q_Ill_Pilon_minMQ_1_minDP_5_Fix_All_Breaks_VCs_All_WiStratificationBy_Pmap_and_SV_4sets_V2:
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







































