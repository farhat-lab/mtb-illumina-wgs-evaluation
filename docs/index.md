This website is related to the repository of for [*Genomic sequence characteristics and the empiric accuracy of short-read sequencing, 2021, BioRxiv*](https://www.biorxiv.org/content/10.1101/2021.04.08.438862v1).


# Results

### Useful Genome-wide statistics and visualizations (across [H37Rv]((https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)), the Mtb reference genome)
From this work we present many useful results that can help guide future genomics studies of the Mtb genome using Illumina WGS. 

[**Full Screen Browser**](jbrowse2/index.html)

<iframe
    src="./jbrowse2/index.html"
    style="border: 1px solid black"
    width="100%"
    height="500"
>
</iframe>




### [Genome Masking Schemes](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/References/Mtb_H37Rv_MaskingSchemes)
1) Refined Low Confidence (RLC) Regions ([RLC_Regions.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.H37Rv.bed))

2) Low Pileup Mappability Regions (K=50bp, E=4 mismatches, [201027_PMap_K50E4_Regions_BELOW_1.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/References/Mtb_H37Rv_MaskingSchemes/201027_PMap_K50E4_Regions_BELOW_1.bed))

3) RLC Regions & Low Pileup Mappability Regions **combined** ([RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed](https://raw.githubusercontent.com/farhat-lab/mtb-illumina-wgs-evaluation/main/Results/B_Extra_UsefulDataFiles/F_Defining_RLC_Regions/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed))
