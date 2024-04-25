#-------------------------------------
# Library
#-------------------------------------

library(Seurat)
library(singleCellTK)
library(data.table)
library(SingleCellExperiment)
library(scuttle)

#-------------------------------------
# Data
#-------------------------------------

SODIR     <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
n         <- length(SODIR)

#-------------------------------------
# Convertion
#-------------------------------------

sce <- readRDS(SODIR)

## give the adequate format for the normalisation
sce <- logNormCounts(sce) 

#-------------------------------------
#Output
#-------------------------------------

exportSCE(sce,
          samplename = "Convert",
          directory  = "03_Input/Atlas",
          type       = "atlas",
          format     = c("Seurat")
          )

