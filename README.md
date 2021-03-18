# RegionalRetina
scripts for regional retina manuscript, which compares the foveal, parafoveal, and macular retinal gene expression at the bulk and single-cell RNA levels

# Overview:

## scRNA-seq analysis
0. For simple questions such as "what cell population expresses me gene?" or "is my gene differentially expressed between the fovea and parafovea?" we point interested users to https://singlecell-eye.org. This web address houses our interactive single-cell RNA sequencing tool Spectacle. Data from this study is fully available at this address.
1. The R script `RegionalRetina_aggregating_filtering_classifying.Rmd` overviews how to take this experiment's count matrices (from GEO), aggregate/normalize/filter the samples, and classify resulting clusters into retinal cell populations. Alternatively, the pre-processed Seurat RData object used for this study is available at   https://drive.google.com/drive/folders/1HWB7T4PqWSyYeBX3uVdAuRCmzYqeDXgA?usp=sharing
2. The R script `github_regional_retina_figures.Rmd` uses the RData object created in (1) to re-create all figures included in this manuscript. The R script `github_plotting_functions.R` is a helper script containing in-house plotting functions for scRNA-seq data. 

## bulk RNA sequencing analysis
The R script `github_bulk_rnaseq_analysis.Rmd` uses edgeR to find differenitally expressed genes between the foveal, parafoveal, and macular retina from the bulk RNA sequencing data. The R script `count_reads.R` overviews how we counted reads after mapping with the aligner STAR. 
