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
library(Seurat)
library(zellkonverter)
library(singleCellTK)

#-------------------------------------
# Data
#-------------------------------------

SCEDIR    <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
n         <- length(SCEDIR)

#-------------------------------------
# Normalisation
#-------------------------------------

for(i in 1:n)
{
  so <- readRDS(file.path(SCEDIR[i]))

  # convert to SingleCellExperiment
  sce_dataset <- as.SingleCellExperiment(so)

  clusts_sce_dataset <- scran::quickCluster(sce_dataset,
                                            method = "igraph",
                                            use.ranks = FALSE, # suggested by the authors
                                            min.size = 100) # require at least 100 cells per cluster

  min.clust = min(table(clusts_sce_dataset))/2
  new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
  sce_dataset = computeSumFactors(sce_dataset, clusters = clusts_sce_dataset, sizes = new_sizes, max.cluster.size = 3000)
    
  sce_dataset <- logNormCounts(sce_dataset, size.factors = sizeFactors(sce_dataset))

  so[["RNA"]] <- SetAssayData(so[["RNA"]],
                              layer = "data", 
                              new.data = logcounts(sce_dataset))

  so$sizeFactors <- sizeFactors(sce_dataset)

  UpdateSeuratObject(so)

  #-------------------------------------
  # Output
  #-------------------------------------
  saveRDS(so, file = file.path(OUTPUTDIR[i]))
  }