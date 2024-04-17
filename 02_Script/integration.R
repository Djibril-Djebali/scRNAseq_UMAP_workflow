#-------------------------------------
# Library
#-------------------------------------

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

#-------------------------------------
# Data
#-------------------------------------

SODIR     <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
n         <- length(SODIR)

#-------------------------------------
# Integration
# Method: https://satijalab.org/seurat/archive/v4.3/integration_mapping
#-------------------------------------

#Create list with our Seurat object to integrate with variable feature selection
so_list <- c()

for(i in 1:n)
{
    so_list    <- append(so_list, readRDS(SODIR[i]))
    so_list[i] <- FindVariableFeatures(so_list[[i]], 
                                       selection.method = "vst", 
                                       nfeatures = 2000,
                                       verbose = FALSE
                                       )
    }

#Rename "celltype" instead of labels the metadata of celltype
for(i in 2:n)
{
    so <- so_list[[i]]
    colnames(so@meta.data)[colnames(so@meta.data) == "labels"] <- "celltype"
    so_list[[i]] <- so
    }

#Identifying anchors,
ANCHORSET <- FindIntegrationAnchors(object.list = so_list[1:n], dims = 1:30)

#Data integration
integrated <- IntegrateData(anchorset = ANCHORSET, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

#-------------------------------------
# Output
#-------------------------------------

saveRDS(integrated, OUTPUTDIR)
