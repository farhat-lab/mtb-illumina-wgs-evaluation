



###### PACBIO ONLY PART ######


rule merge_Wi_CAT_All_PacBio_FQ_GZs:
    input:
        lambda wildcards: expand("{i_PacBio_FQ_PATH}", i_PacBio_FQ_PATH = SampleID_To_PacBio_FQs_Dict[wildcards.sampleID_WiPacBio])
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    threads: 1
    shell:
        "cat {input} > {output}"




###################################################
############# NanoPlot (Long read QC) #############
###################################################

rule nanoplot_QC:
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Nanoplot_QC/{sampleID_WiPacBio}.subreads.NanoPlot/NanoPlot-report.html"
    conda:
        "CondaEnvs/nanoplot_128_Conda.yml"
    threads: 4
    params:
        NanoPlot_OutputDir_PATH = output_Dir + "/{sampleID_WiPacBio}/pacbio/Nanoplot_QC/{sampleID_WiPacBio}.subreads.NanoPlot/"
    shell:
        "NanoPlot -t {threads} --fastq {input} -o {params.NanoPlot_OutputDir_PATH}"


##########################################################
############# TSV of PacBio read lengths ###############
##########################################################


rule PacBio_Subreads_GetReadLengthsTSV:
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.ReadLengths.tsv"
    params:
        bioawk_GetReadLengths_TEXT = "'{print $name, length($seq)}'"
    threads: 1
    shell:
        "bioawk -c fastx {params.bioawk_GetReadLengths_TEXT} {input} > {output}"











rule flye_Assemble: # Flye v2.6 w/ asmCov = 200
    input:
        pb_subreads_fq = output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    output:
        assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta",
        assembly_info_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly_info.txt"
    conda:
       "CondaEnvs/PacBio_Software_py27_Conda.yml"
    threads: 10
    params:
        Flye_OutputDir_PATH = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/"
    shell:
        "flye --pacbio-raw {input.pb_subreads_fq} --out-dir {params.Flye_OutputDir_PATH} "
        "--genome-size 5m --threads {threads} --asm-coverage 200 --iterations 3"




########################################################################
######### Generate list of circular and non-circular contigs ###########
########################################################################
# https://stackoverflow.com/questions/50965417/awk-command-fails-in-snakemake-use-singularity

#rule output_CircularContigs:
#    input:
#        assembly_info_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly_info.txt"
#    output:
#        assembly_circularcontigs_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/circularContigs_List.txt",
#        assembly_not_circularcontigs_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/NoncircularContigs_List.txt"
#    threads: 1
#    shell:
#        "awk '$4 == \"Y\" {{print $1}}' {input.assembly_info_txt} > {output.assembly_circularcontigs_txt} \n"
#        "awk '$4 == \"N\" {{print $1}}' {input.assembly_info_txt} > {output.assembly_not_circularcontigs_txt}"

rule output_CircularContigs:
    input:
        assembly_info_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly_info.txt"
    output:
        assembly_circularcontigs_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/circularContigs_List.txt",
        assembly_not_circularcontigs_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/NoncircularContigs_List.txt"
    threads: 1
    shell:
        "awk '$4 == \"+\" {{print $1}}' {input.assembly_info_txt} > {output.assembly_circularcontigs_txt} \n"
        "awk '$4 == \"-\" {{print $1}}' {input.assembly_info_txt} > {output.assembly_not_circularcontigs_txt}"





###################################################################################
######### CIRCLATOR for setting start at DnaA (Assuming Circular genome) ##########
###################################################################################

rule circlator_FixStart_DnaA:
    input:
        #assembly_not_circularcontigs_txt = output_Dir + "{sampleID}/pacbio/Flye_Assembly/NoncircularContigs_List.txt",
        DnaA_Seq_fa = H37rv_DnaA_FA_PATH,
        flye_assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta"
    output:
        flye_assembly_FixStart_assembly = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    conda:
        "CondaEnvs/circulator_151py37_Conda.yml"
    threads: 1
    shell:
        "circlator fixstart --genes_fa {input.DnaA_Seq_fa} {input.flye_assembly_fa} {output_Dir}/{wildcards.sampleID_WiPacBio}/pacbio/Flye_Assembly/{wildcards.sampleID_WiPacBio}.flyeassembly.fixstart"





rule samtools_faidx_FlyeAssembly:
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta.fai"
    conda:
        "CondaEnvs/PacBio_Software_py27_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"




rule CP_PacBio_FlyeAssembly_To_I3_Dir: # This is a TEMP fix to just skip the quiver polishing steps
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.fixstart.fasta",
    shell: 
        "cp {input} {output}"

### Fitler PacBio assembly (GC3) for only contigs greater than 100kb

### To Do:
### 1) Filter out small contigs (< 100kb) from the PacBio.GC3 (PacBio polished assembly)
### 2) Rerun all analysis using only the large/completed genome assembly

rule filterByLength_100kbContigs_FlyeAssembly_I3:
    input:
        PacBio_Flye_Assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.fixstart.fasta",
    output:
        PacBio_Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta",
        PacBio_Flye_Assembly_Renamed_100KbContigs_FAI = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta.fai"
    threads: 1
    conda:
        "CondaEnvs/bioinfo_util_env_V1.yml" # "CondaEnvs/bioinfo_util_env.yml"
    shell:
        " bioawk -c fastx '{{ print \">{wildcards.sampleID_WiPacBio}_\"$name \"\\n\" $seq }}' {input.PacBio_Flye_Assembly_fa} "
        " | "
        " bioawk -c fastx '{{ if(length($seq) > 100000) {{ print \">\"$name; print $seq }}}}' > {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} \n"
        " samtools faidx {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} "



rule FastANI_FlyeAssembly_I3:
    input:
        PacBio_Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta",
        H37rv_Ref_fa = refGenome_FA_PATH,     
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/FastANI_OutputDirs/FlyeAssembly_I3_PBonly_FastANI/{sampleID_WiPacBio}.flyeassembly.I3.AssemblyToH37rv.FastANI.txt"
    conda:
        "CondaEnvs/fastani_1_3_0_Conda.yml"
    shell:
        "fastANI -q {input.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} -r {input.H37rv_Ref_fa} -o {output}"



###########################################
####### Prokka Genome Annotation ##########
###########################################


rule Prokka_Anno_FlyeAssembly_I3_PB_DeNovo:
    input:
        PacBio_Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta",
    output:
        Q3_Prokka_Anno_FNA = output_Dir + "/{sampleID_WiPacBio}/Prokka/Flye_I3_PB_DeNovo_Prokka_Anno/{sampleID_WiPacBio}.fna"
    conda:
        "CondaEnvs/Prokka_1130_Conda.yml" #"CondaEnvs/Prokka_1145_Conda.yml"
    threads: 8
    shell: 
        "prokka --force --cpus {threads} --genus Mycobacterium --species tuberculosis --strain {wildcards.sampleID_WiPacBio}"
        " --addgenes --mincontiglen 100000 --centre XXX --proteins {H37rv_GBK_PATH} "
        " --outdir {output_Dir}/{wildcards.sampleID_WiPacBio}/Prokka/Flye_I3_PB_DeNovo_Prokka_Anno/ --prefix {wildcards.sampleID_WiPacBio} {input.PacBio_Flye_Assembly_Renamed_100KbContigs_FA}"






########################################################
### LINEAGE CALLING script using VCF as input ###
########################################################


rule LineageCall_Coll_Flye_I3_PP_AlignTo_H37rv:
    input:
        MM2_Flye_I3_PP_To_H37rv_VCF = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PP_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.mm2.Flye_I3_PP_AssemblyToH37rv.paftools.vcf",
        Coll_SNPs_PATH = "References/LineageCalling_SNP_Schemes/coll.tsv"
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/LineageCalling/LineageCall_Flye_I3_PP_MM2_AlignTo_H37rv/{sampleID_Wi_Ill_And_PB}.AtoRef.FlyeI3_PP.lineage_call.Coll.txt"
    shell:
        "Scripts/fast-lineage-caller-vcf.py {input.MM2_Flye_I3_PP_To_H37rv_VCF} {input.Coll_SNPs_PATH} > {output} "




rule LineageCall_Coll_GC3_PBonly_AlignTo_H37rv:
    input:
        MM2_Flye_I3_PBonly_To_H37rv_VCF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Minimap2_Flye_I3_PBonly_AlignTo_H37rv/{sampleID_WiPacBio}.mm2.Flye_I3_PBonly_AssemblyToH37rv.vcf",
        Coll_SNPs_PATH = "References/LineageCalling_SNP_Schemes/coll.tsv"
    output:
        output_Dir + "/{sampleID_WiPacBio}/LineageCalling/LineageCall_Flye_I3_PBonly_MM2_AlignTo_H37rv/{sampleID_WiPacBio}.AtoRef.Flye_I3_PBonly.lineage_call.Coll.txt"
    shell:
        "Scripts/fast-lineage-caller-vcf.py {input.MM2_Flye_I3_PBonly_To_H37rv_VCF} {input.Coll_SNPs_PATH} > {output} "














#######################################################################
######## Align PacBio subreads to H37Rv reference w/ Minimap2 #########
#######################################################################

rule align_PacBio_Subreads_To_H37Rv_Minimap2:
    input:
        pb_subreads_fq = output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
        Ref_Genome_FA = refGenome_FA_PATH,
    output:
        pb_subreads_MM2_To_H37Rv_SAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.sam"
    threads: 8
    shell:
        "minimap2 -t {threads} -ax map-pb {input.Ref_Genome_FA} {input.pb_subreads_fq} > {output.pb_subreads_MM2_To_H37Rv_SAM}"
        # "minimap2 --MD -t {threads} -ax map-pb {input.Ref_Genome_FA} {input.pb_subreads_fq} > {output.pb_subreads_MM2_To_H37Rv_SAM}"



rule samtools_ViewAndSort_PacBio_Subreads_To_H37Rv_Minimap2:
    input:
        pb_subreads_MM2_To_H37Rv_SAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.sam",
    output:
        pb_subreads_MM2_To_H37Rv_BAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam",
        pb_subreads_MM2_To_H37Rv_BAM_BAI = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.bai",
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input.pb_subreads_MM2_To_H37Rv_SAM} | samtools sort - > {output.pb_subreads_MM2_To_H37Rv_BAM} \n"
        "samtools index {output.pb_subreads_MM2_To_H37Rv_BAM}"




rule samtools_Depth_PacBio_Subreads_To_H37Rv_Minimap2:
    input:
        pb_subreads_MM2_To_H37Rv_BAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam",
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.depth.txt"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools depth -a {input} > {output}"



rule samtools_Depth_AverageAll_PacBio_Subreads_To_H37Rv_Minimap2:
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.depth.txt"
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.depth.averaged.txt"
    shell:
        "awk '{{sum+=$3}} END {{ print \"Average = \",sum/NR}}' {input} > {output}"




rule pilon_VarCalling_PacBio_Subreads_AlignTo_H37rv_DefaultParam:
    input:
        Ref_fa = refGenome_FA_PATH,
        pb_subreads_MM2_To_H37Rv_BAM = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam",
        pb_subreads_MM2_To_H37Rv_BAM_BAI = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/PacBio_Subreads_AlignedToH37Rv_Minimap2/{sampleID_WiPacBio}.pb.subreads.H37Rv.minimap2.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Pilon_PacBio_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiPacBio}.PacBio.AlnToH37rv.vcf",
        pilon_ChangesFile = output_Dir+ "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Pilon_PacBio_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiPacBio}.PacBio.AlnToH37rv.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Pilon_PacBio_AlignedTo_H37rv_DefaultParameters_VariantCalling/"
    conda:
        "CondaEnvs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx85g --genome {input.Ref_fa} --pacbio  {input.pb_subreads_MM2_To_H37Rv_BAM} --output {wildcards.sampleID_WiPacBio}.PacBio.AlnToH37rv"
        " --outdir {params.Pilon_OutputDir_PATH} --changes --tracks --fix snps,indels --vcf"




rule calculate_F2_Score_PacBio_Subreads:
    input:
        pilon_VCF = output_Dir+ "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/Pilon_PacBio_AlignedTo_H37rv_DefaultParameters_VariantCalling/{sampleID_WiPacBio}.PacBio.AlnToH37rv.vcf",
        lineage_def_ref_pos_TXT = "Scripts/positions_to_select.txt",
        Coll2014_LinSpeSNPs_final_CSV = "Scripts/Coll2014_LinSpeSNPs_final.csv",        
    output:
        F2_Score_TXT = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/F2_Calculation/{sampleID_WiPacBio}_F2.txt",
        #TEMP_BCF = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf",
        #TEMP_BCF_CSI = output_Dir + "/{sampleID_WiIll}/IlluminaWGS/F2_Calculation/current.bcf",
    params:
        F2_OutputDir_PATH = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/F2_Calculation",
        TEMP_BCF = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/F2_Calculation/current.bcf",
        TEMP_BCF_CSI = output_Dir + "/{sampleID_WiPacBio}/pacbio_VariantCallingVersusH37Rv/F2_Calculation/current.bcf.csi",
    shell:
        "python Scripts/lineage_defining_SNP_depths_collection_and_F2_calculation_PacBio.py {input.pilon_VCF} {params.F2_OutputDir_PATH} {input.lineage_def_ref_pos_TXT} {wildcards.sampleID_WiPacBio} {input.Coll2014_LinSpeSNPs_final_CSV} \n"
        "rm {params.TEMP_BCF} {params.TEMP_BCF_CSI} "





















### DONE W/ PACBIO ONLY PART OF ANALYSIS ###





