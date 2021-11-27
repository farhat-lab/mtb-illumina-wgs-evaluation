# Evaluation of short-read sequencing performance across the Mtb genome<a name="evaluation-of-short-read-sequencing-performance-across-the-mtb-genome"></a>

This repository contains the bioinformatics pipeline and code needed to
reproduce the analysis for
[*Genomic sequence characteristics and the empiric accuracy of short-read sequencing, 2021, BioRxiv*](https://www.biorxiv.org/content/10.1101/2021.04.08.438862v1).
All necessary bioinformatics processing is implemented using the
[SnakeMake](https://snakemake.github.io/) workflow system, along with all
downstream analysis available in Jupyter notebooks using Python 3.7. All
software dependencies are defined for processing and analysis steps using
[Conda](https://docs.conda.io/en/latest/).

Go [here](https://farhat-lab.github.io/mtb-illumina-wgs-evaluation/jbrowse2/index.html)
to interactively explore the key results from this paper in an online genome
browser.

The provided Snakemake pipeline implements all steps related to Assembly,
alignment, and other key processing steps.

## Table of Contents<a name="table-of-contents"></a>

<!-- mdformat-toc start --slug=github --maxlevel=6 --minlevel=1 -->

- [Evaluation of short-read sequencing performance across the Mtb genome](#evaluation-of-short-read-sequencing-performance-across-the-mtb-genome)
  - [Pre-requisites](#pre-requisites)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Supporting Data Analysis](#supporting-data-analysis)
  - [Results](#results)
    - [Useful genome-wide statistics and visualizations for H37Rv (e.g., the Mtb reference genome](#useful-genome-wide-statistics-and-visualizations-for-h37rv-eg-the-mtb-reference-genome)
    - [Pileup Mappability (K = 50 bp, E >= 4 mismatches)](#pileup-mappability-k--50-bp-e--4-mismatches)
    - [Pileup Mappability (K = 100 bp, E >= 4 mismatches)](#pileup-mappability-k--100-bp-e--4-mismatches)
    - [Empirical Base pair Recall (EBR, 36 clinical Mtb isolates)](#empirical-base-pair-recall-ebr-36-clinical-mtb-isolates)
  - [Interactive EBR & Pileup Mappability visualization](#interactive-ebr--pileup-mappability-visualization)
  - [Genome Masking Schemes](#genome-masking-schemes)
  - [License](#license)

<!-- mdformat-toc end -->

## Pre-requisites<a name="pre-requisites"></a>

- A TSV with sample metadata and input FASTQ paths (Illumina and PacBio).

Find examples [here](./Data/201202_PMP_SM_50CI_AllDataSets_InputSeqDataPaths) directory.

- Path to a target output directory (e.g., `/home/user/snakemake-output`)

- A [config file](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/blob/main/Snakemake_Rules_And_Config/config_PMP_V6.json)
  with file path to H37Rv reference files
  [`NC_000962.3`](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3) in FASTA and
  GBK formats.

- A [cluster config file](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/blob/main/Snakemake_Rules_And_Config/clusterConfig_PMP_V10.json)
  specifying computational resource requirements for steps of the pipeline when
  using [SLURM](https://slurm.schedmd.com/documentation.html)

## Installation<a name="installation"></a>

All dependencies needed to reproduce this analysis can be installed via
[Conda](https://docs.conda.io/en/latest/).

1. Clone the repository

```bash
git clone https://github.com/farhat-lab/mtb-illumina-wgs-evaluation
cd mtb-illumina-wgs-evaluation/
```

2. Create a conda environment named `CoreEnv_MtbEval_V1`

```bash
conda env create --file CondaEnvs/Core_CondaEnv_Mtb_Illumina_WGS_Benchmarking_V1.yml -n CoreEnv_MtbEval_V1
```

This command creates a Conda environment containing the following:

- Snakemake workflow engine
- Bioinformatics software dependencies
- Python 3 libraries for downstream analysis

3. Activate environment

```bash
conda activate CoreEnv_MtbEval_V1
```

## Usage<a name="usage"></a>

To run the SnakeMake workflow, run the following

1. Define path to the target output directory

```bash
OUTPUT_DIR="../Mtb_PacBio_And_Illumina_Manuscript_SnakeMake_Output_V1"
```

2. Define config file

```bash
SNAMEMAKE_CONFIG_FILE="Snakemake_Rules_And_Config/config_PMP_V6.json"
```

3. Define SLURM config file

```bash
SLURM_CONFIG_FILE="Snakemake_Rules_And_Config/clusterConfig_PMP_V10.json"
```

4. Define a path to TSV that specifies Sample Metadata and FASTQ PATHs

```bash
INPUT_DATA_TSV_DIR="Data/201202_PMP_SM_50CI_AllDataSets_InputSeqDataPaths"
INPUT_SAMPLEINFO_TSV="${INPUT_DATA_TSV_DIR}/201202_MTb_50CI_Peru_ChinerOms_Ngabonziza_TBPortals_PacBioDatasetsMerged_SampleInfo_InputFQs.tsv"
```

5. Create logs directory

```bash
mkdir -p "${OUTPUT_DIR}/logs/cluster"
```

6. Submit job to Slurm to run the SnakeMake pipeline

```bash
--cluster 'sbatch \
-p {cluster.p} \
-n {cluster.n} \
-t {cluster.t} \
--mem {cluster.mem} \
-o "${OUTPUT_DIR}/{cluster.o}" \
-e "${OUTPUT_DIR}/{cluster.e}"'
```

## Supporting Data Analysis<a name="supporting-data-analysis"></a>

The
[DataAnalysis/](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/DataAnalysis)
directory contains Jupyter notebooks for downstream data processing as well as table
and figure generation related to this work.

## Results<a name="results"></a>

### Useful genome-wide statistics and visualizations for [H37Rv](<(https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)>) (e.g., the Mtb reference genome<a name="useful-genome-wide-statistics-and-visualizations-for-h37rv-eg-the-mtb-reference-genome"></a>

We present many valuable results from this work that can help guide future
genomics studies of the Mtb genome using Illumina WGS.

### Pileup Mappability (K = 50 bp, E >= 4 mismatches)<a name="pileup-mappability-k--50-bp-e--4-mismatches"></a>

[201027_H37rv_PileupMappability_K50_E4.bedgraph](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/201027_H37rv_PileupMappability_K50_E4.bedgraph)

[201027_H37rv_PileupMappability_K50_E4.bw](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/201027_H37rv_PileupMappability_K50_E4.bw)

### Pileup Mappability (K = 100 bp, E >= 4 mismatches)<a name="pileup-mappability-k--100-bp-e--4-mismatches"></a>

[201027_H37rv_PileupMappability_K100_E4.bedgraph](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/201027_H37rv_PileupMappability_K100_E4.bedgraph)

[201027_H37rv_PileupMappability_K100_E4.bw](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/201027_H37rv_PileupMappability_K100_E4.bw)

### Empirical Base pair Recall (EBR, 36 clinical Mtb isolates)<a name="empirical-base-pair-recall-ebr-36-clinical-mtb-isolates"></a>

[210112_EBR_V7_36CI.bedgraph](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/210112_EBR_V7_36CI.bedgraph)

[210112_EBR_V7_36CI.bw](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/C_BrowserTracks/210112_EBR_V7_36CI.bw)

## [Interactive EBR & Pileup Mappability visualization](https://farhat-lab.github.io/mtb-illumina-wgs-evaluation/jbrowse2/index.html)<a name="interactive-ebr--pileup-mappability-visualization"></a>

These results can easily be explored in a browser-based JBroswe2 genome browser,
found
[here](https://farhat-lab.github.io/mtb-illumina-wgs-evaluation/jbrowse2/index.html).
If you want to visualize the genome-wide statistics in your genome browser,
use the `.bedgraph` and `.bw` files listed above.

## [Genome Masking Schemes](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/References/Mtb_H37Rv_MaskingSchemes)<a name="genome-masking-schemes"></a>

1. Refined Low Confidence (RLC) Regions
   ([RLC_Regions.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.H37Rv.bed))

1. Low Pileup Mappability Regions (K=50bp, E=4 mismatches,
   [201027_PMap_K50E4_Regions_BELOW_1.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed))

1. RLC Regions & Low Pileup Mappability Regions **combined**
   ([RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/B_Extra_UsefulDataFiles/F_Defining_RLC_Regions/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed))

## License<a name="license"></a>

This repository is distributed under the [MIT license terms](LICENSE).
