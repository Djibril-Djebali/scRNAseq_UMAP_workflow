#-------------------------------------
# Library
#-------------------------------------
library(singleCellTK)
library(data.table)
library(SingleCellExperiment)

#-------------------------------------
# Data
#-------------------------------------

SODIR     <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
n         <- length(SODIR)

#-------------------------------------
# Convertion
#-------------------------------------
for(i in 1:n){
    so  <- readRDS(SODIR[i])
    sce <- convertSeuratToSCE(so)

    #-------------------------------------
    #Output
    #-------------------------------------

    saveRDS(sce, file = file.path(OUTPUTDIR[i]))
    }

