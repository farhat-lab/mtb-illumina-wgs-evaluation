# Evaluation of short-read sequencing performance across the Mtb genome

This repository contains the bioinformatics pipeline and code needed to reproduce the analysis for [*Genomic sequence characteristics and the empiric accuracy of short-read sequencing, 2021, BioRxiv*](https://www.biorxiv.org/content/10.1101/2021.04.08.438862v1). All key bioinformatics processing is implemented using the [SnakeMake](https://snakemake.github.io/) workflow language, along with all downstream analysis available in Jupyter notebooks using Python 3.7. All software dependencies are defined for processing and analysis steps using [Conda](https://docs.conda.io/en/latest/). 


## Contents
- [Installation](#Installation)
- [SnakeMake Pipeline]()
- [Data Analysis]()
- [License](#License)


## A) Installation
All dependencies needed to reproduce analysis can be installed via Conda(https://docs.conda.io/en/latest/) .
```
# 1) Clone repository
git clone https://github.com/farhat-lab/mtb-illumina-wgs-evaluation

# 2) Create a conda environment named 'CoreEnv_MtbEval_V1'
cd mtb-illumina-wgs-evaluation/

conda env create --file envname.yml -n CoreEnv_MtbEval_V1

# 3) Activate environment (used for SnakeMake pipeline and data analysis)
conda activate CoreEnv_MtbEval_V1
```

The above command will create a Conda environment with:
- a) the Snakemake workflow engine,
- b) bioinformatics software dependencies, 
- c) Python libraries for downstream analysis.


## B) Running the SnakeMake Pipeline (Processing of Illumina and PacBio data)

The provided Snakemake pipeline implements all steps related to Assembly, alignment, and other key processing steps.

To run the SnakeMake pipeline you need to provide:
- a) Path to a target outout directory
- b) A config file with paths to H37Rv reference files (FASTA and GBK formats) 
- c) A config file specifing computational resource requirements for 


``` 
# Define path to the target output directory
targetOutput_Dir="../Mtb_PacBio_And_Illumina_Manuscript_SnakeMake_Output_V1"

SnakeMake_ConfigFile="Snakemake_Rules_And_Config/config_PMP_V6.json"

SLURM_Cluster_Config="Snakemake_Rules_And_Config/clusterConfig_PMP_V10.json"

mkdir ${targetOutput_Dir}

# Define path to TSV that specifies Sample Metadata and FASTQ PATHs

inputData_TSV_Dir="Data/201202_PMP_SM_50CI_AllDataSets_InputSeqDataPaths"

input_SampleInfo_TSV="${inputData_TSV_Dir}/201202_MTb_50CI_Peru_ChinerOms_Ngabonziza_TBPortals_PacBioDatasetsMerged_SampleInfo_InputFQs.tsv"


mkdir -p ${targetOutput_Dir}/O2logs/cluster/

snakemake -s SnakeFile_Main_Processing.smk.py --config output_dir=${targetOutput_Dir} inputSampleData_TSV=${input_SampleInfo_TSV} --configfile ${inputConfigFile} -p --use-conda -j 50 --cluster-config  ${SLURM_Cluster_Config}  --cluster "sbatch -p {cluster.p} -n {cluster.n}  -t {cluster.t} --mem {cluster.mem} -o ${targetOutput_Dir}/{cluster.o} -e ${targetOutput_Dir}/{cluster.e}" --latency-wait 35 -k 
``` 

## C) Supporting Data Analysis 

The [DataAnalysis/](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/DataAnalysis) directory contains Jupyter notebooks for downstream data processing, table generation, and figure generation.





## D) Useful results for future Mtb genomics analyses

### Useful Genome-wide statistics and visualizations (H37Rv, the Mtb reference genome)
From this work we present many useful results that can help guide future genomics studies of the Mtb genome using Illumina WGS. 


### Genome Masking Schemes
### PROVIDED GENOME MASKING SCHEMES HERE




## License
This repository is distributed under the [MIT license terms](LICENSE).

