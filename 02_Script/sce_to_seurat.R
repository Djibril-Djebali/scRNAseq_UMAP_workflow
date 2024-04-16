#-------------------------------------
# Library
#-------------------------------------
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
#but no idea if it will cause cause an issue in next step?
sce <- logNormCounts(sce) 

so  <- convertSCEToSeurat(sce)

#-------------------------------------
#Output
#-------------------------------------

saveRDS(so, file = file.path(OUTPUTDIR))