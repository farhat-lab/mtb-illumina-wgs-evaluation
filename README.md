# Evaluation of short-read sequencing performance across the Mtb genome

This repository contains the bioinformatics pipeline and code needed to reproduce the analysis for [*Genomic sequence characteristics and the empiric accuracy of short-read sequencing, 2021, BioRxiv*](https://www.biorxiv.org/content/10.1101/2021.04.08.438862v1). All key bioinformatics processing is implemented using the [SnakeMake](https://snakemake.github.io/) workflow language, along with all downstream analysis available in Jupyter notebooks using Python 3.7. All software dependencies are defined for processing and analysis steps using [Conda](https://docs.conda.io/en/latest/). 


## Contents
- [Installation](#Installation)
- [Key Results](#Installation)
- [Useful References and Visualizations]()
- [SnakeMake Pipeline]()
- [Data Analysis]()
- [License](#License)

## Key Results 

In this work we used long read sequencing to throughly benchmark Illumina sequencing performance across the Mtb genome. From this work we produced many useful insights for future studies of the Mtb genome using Illumina WGS.

## Useful Genome-wide statistics and visualizations (H37Rv, the Mtb reference genome)
From this work we present many useful results that can help guide future genomics studies of the Mtb genome using Illumina WGS. 


## Genome Masking Schemes
### PROVIDED GENOME MASKING SCHEMES HERE



## Installation
All dependencies needed to reproduce analysis can be installed via Conda(https://docs.conda.io/en/latest/) .
```
```
The above command will create a Conda environment with:
- a) the Snakemake workflow engine,
- b) bioinformatics software dependencies, 
- c) Python libraries for downstream analysis.



## A) Bioinformatics Pipeline (Processing Illumina and PacBio data)

The Snakemake pipeline implements all steps related to Assembly, alignment, and other key processing steps.
```
CODE HERE FOR RUNNING SNAKEMAKE WORKFLOW
```

## B) Supporting Data Analysis 

The ___ directory contains Jupyter notebooks for downstream data processing, table generation, and figure generation.
The direcrory structure is as follows.



## License
This repository is distributed under the [MIT license terms](LICENSE).

