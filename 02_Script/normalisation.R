#-------------------------------------
# Normalisation : we follow the process from Pijuan-Sala et al. 2019 
# (doi:10.1038/s41586-019-0933-9)
#-------------------------------------

#-------------------------------------
# Library
#-------------------------------------
library(data.table)
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(Matrix)
library(irlba)
library(zellkonverter)
library(singleCellTK)

#-------------------------------------
# Data
#-------------------------------------
INPUTDIR = snakemake@input[["input"]]

DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "03_Output")

SAMPLE_ID = snakemake@params[["sample_list"]]
MOUSEGASTRULATION_SAMPLES =  snakemake@params[["mousegastrulation_samples"]]

#-------------------------------------
# Normalisation
#-------------------------------------
sce <- readRDS(INPUTDIR)

## STEP from https://rstudio-pubs-static.s3.amazonaws.com/699579_c6be4bf3220746088bdfd12a61aa15c4.html
# and https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/4_normalisation/normalise.Rmd
# For pre-clustering, we use scran's `quickCluster` function, using the `igraph` method. We specify a maximum cluster size of 3000 cells and a minimum cluster size of 100 cells.
clusts <- as.numeric(quickCluster(sce,
                        method = "igraph",
                        use.ranks = FALSE, # suggested by the authors
                        min.size = 100)) # require at least 100 cells per cluster
# Number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

## Use scuttle to normalize function, which calculates log2-transformed normalized expression values.
# This is done by dividing each count by its size factor, adding a pseudo-count and log-transforming.
sce <- logNormCounts(sce)

#-------------------------------------
# Output
#-------------------------------------
#saveRDS(sce, file = file.path("03_Output/Normalised/mouse_gastrulation_normdata.rds"))

exportSCE(
  sce,
  samplename = "Normalised",
  type       = "mouse_gastrulation",
  directory  = "03_Output",
  format     = c("Seurat")
  )