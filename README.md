# Evaluation of short-read sequencing performance across the Mtb genome

This repository contains the bioinformatics pipeline and code needed to reproduce the analysis for [*Genomic sequence characteristics and the empiric accuracy of short-read sequencing, 2021, BioRxiv*](https://www.biorxiv.org/content/10.1101/2021.04.08.438862v1). All key bioinformatics processing is implemented using the [SnakeMake](https://snakemake.github.io/) workflow system, along with all downstream analysis available in Jupyter notebooks using Python 3.7. All software dependencies are defined for processing and analysis steps using [Conda](https://docs.conda.io/en/latest/). 


## Contents
- [Installation](#Installation)
- [SnakeMake Pipeline](#Processing-of-Illumina-and-PacBio-sequencing)
- [Data Analysis](#Supporting-Data-Analysis)
- [Results](#Results)
- [License](#License)


# Installation
All dependencies needed to reproduce this analysis can be installed via [Conda](https://docs.conda.io/en/latest/) .
```
# 1) Clone repository
git clone https://github.com/farhat-lab/mtb-illumina-wgs-evaluation

# 2) Create a conda environment named 'CoreEnv_MtbEval_V1'
cd mtb-illumina-wgs-evaluation/

conda env create --file CondaEnvs/Core_CondaEnv_Mtb_Illumina_WGS_Benchmarking_V1.yml -n CoreEnv_MtbEval_V1

# 3) Activate environment (used for SnakeMake pipeline and data analysis)
conda activate CoreEnv_MtbEval_V1
```

The above command will create a Conda environment with:
- a) the Snakemake workflow engine,
- b) bioinformatics software dependencies, 
- c) Python libraries for downstream analysis.


# Processing of Illumina and PacBio sequencing 

The provided Snakemake pipeline implements all steps related to Assembly, alignment, and other key processing steps.

To run the SnakeMake pipeline you need to provide:
- a) A TSV with sample metadata and input FASTQ paths (Illumina and PacBio). Examples can be found in `Data/201202_PMP_SM_50CI_AllDataSets_InputSeqDataPaths`.
- b) Path to a target outout directory
- c) A [config file](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/blob/main/Snakemake_Rules_And_Config/config_PMP_V6.json) with paths to H37Rv reference files ([NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)) in FASTA and GBK formats.
- d) A [cluster config file](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/blob/main/Snakemake_Rules_And_Config/clusterConfig_PMP_V10.json) specifing computational resource requirements for steps of the pipeline when using [SLURM](https://slurm.schedmd.com/documentation.html)


To run the SnakeMake workflow, run the following bash commands:
``` 
# Define path to the target output directory
targetOutput_Dir="../Mtb_PacBio_And_Illumina_Manuscript_SnakeMake_Output_V1"

# Define config file 
SnakeMake_ConfigFile="Snakemake_Rules_And_Config/config_PMP_V6.json"

# Define SLURM config file
SLURM_Cluster_Config="Snakemake_Rules_And_Config/clusterConfig_PMP_V10.json"

# Define path to TSV that specifies Sample Metadata and FASTQ PATHs

inputData_TSV_Dir="Data/201202_PMP_SM_50CI_AllDataSets_InputSeqDataPaths"
input_SampleInfo_TSV="${inputData_TSV_Dir}/201202_MTb_50CI_Peru_ChinerOms_Ngabonziza_TBPortals_PacBioDatasetsMerged_SampleInfo_InputFQs.tsv"

# Run SnakeMake pipeline
mkdir ${targetOutput_Dir}
mkdir -p ${targetOutput_Dir}/O2logs/cluster/

snakemake -s SnakeFile_Main_Processing.smk.py --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} --configfile ${inputConfigFile} -p --use-conda -j 50 --cluster-config  ${SLURM_Cluster_Config}  --cluster "sbatch -p {cluster.p} -n {cluster.n}  -t {cluster.t} --mem {cluster.mem} -o ${targetOutput_Dir}/{cluster.o} -e ${targetOutput_Dir}/{cluster.e}" --latency-wait 35 -k 
``` 
The final command will begin the SnakeMake pipeline using the SLURM workload manager. 



# Supporting Data Analysis 

The [DataAnalysis/](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/DataAnalysis) directory contains Jupyter notebooks for downstream data processing, table generation, and figure generation related to this work.


# Results

### Useful Genome-wide statistics and visualizations (across [H37Rv]((https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)), the Mtb reference genome)
From this work we present many useful results that can help guide future genomics studies of the Mtb genome using Illumina WGS. 

### [Interactive EBR & Pileup Mappability visualizations](https://farhat-lab.github.io/mtb-illumina-wgs-evaluation/jbrowse2/index.html)



### [Genome Masking Schemes](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/References/Mtb_H37Rv_MaskingSchemes)
1) Refined Low Confidence (RLC) Regions ([RLC_Regions.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.H37Rv.bed))

2) Low Pileup Mappability Regions (K=50bp, E=4 mismatches, [201027_PMap_K50E4_Regions_BELOW_1.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed))

3) RLC Regions & Low Pileup Mappability Regions **combined** ([RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/B_Extra_UsefulDataFiles/F_Defining_RLC_Regions/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed))


## License
This repository is distributed under the [MIT license terms](LICENSE).

