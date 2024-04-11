#-------------------------------------
# Library
#-------------------------------------
library(singleCellTK)
library(data.table)
library(SingleCellExperiment)

#-------------------------------------
# Convertion
#-------------------------------------

SODIR     <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
n         <- length(SODIR)

for(i in 1:n){
    so  <- readRDS(SODIR[i])
    sce <- convertSeuratToSCE(so)

    #-------------------------------------
        #Output
    #-------------------------------------

    saveRDS(so, file = file.path(OUTPUTDIR[i]))
    }

