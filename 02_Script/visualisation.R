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
INPUTDIR <- snakemake@input[["input"]]
OUTPUTDIR <- snakemake@output[["output"]]
integrated <- readRDS(INPUTDIR)

p1 <- DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
      NoLegend()

pdf(OUTPUTDIR, height = 10, width = 20)
p1 + p2
dev.off()