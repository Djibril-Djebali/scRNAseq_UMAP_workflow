#..............................................................
# Mouse gastrulation reference :
# doi:10.1038/s41586-019-0933-9
#..............................................................

#-------------------------------------
# Library
#-------------------------------------

library(data.table)
library(SingleCellExperiment)
BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)
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

SCEDIR    <- snakemake@output[["output"]]

##  Loading the samples from the atlas 
MOUSEGASTRULATION_SAMPLES <- snakemake@params[["mousegastrulation_samples"]]
MOUSEGASTRULATION_SAMPLES <- type.convert(MOUSEGASTRULATION_SAMPLES, dec=".", as.is = TRUE)

# Loading all the samples in a sce object
sce <- EmbryoAtlasData(samples = c(MOUSEGASTRULATION_SAMPLES))

# Change ENSEMBL name of the count matrix by corresponding gene name 
sce <- setRowNames(sce, "SYMBOL")

# Remove the cell with Na as labeltype
sce <- sce[,!is.na(colData(sce)$celltype)]

#-------------------------------------
# Output
#-------------------------------------

saveRDS(sce, file = file.path(SCEDIR))