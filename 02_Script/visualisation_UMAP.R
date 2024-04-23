#===============================================================================
#Library
#===============================================================================

library(ggplot2)
library(ggnewscale)
library(Seurat)
library(cowplot)

#===============================================================================
#Data
#===============================================================================

INPUTDIR   <- snakemake@input[["input"]]
OUTPUTDIR1 <- snakemake@output[["output1"]]
OUTPUTDIR2 <- snakemake@output[["output2"]]
OUTPUTDIR3 <- snakemake@output[["output3"]]
SAMPLEID   <- snakemake@params[["params"]]
integrated <- readRDS(INPUTDIR)

#Create a data frame 4 x number of cells, with each cells assigned the 
#parameters UMAP_1, UMAP_2, origin, celltype
data          <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
data$origin   <- unlist(as.list(integrated@meta.data[["orig.ident"]]))
data$celltype <- unlist(as.list(integrated@meta.data[["celltype"]]))
data$sample   <- unlist(as.list(integrated@meta.data[["sample"]]))
integrated    <- NA

#Split for each conditionsnakemake@output[["output"]]
data_split    <- split(data, data$origin)
n             <- length(data_split)

data_ref  <- data_split[[SAMPLEID[1]]] 
data_allh <- data_split[[SAMPLEID[2]]]  
for(i in 3:(n))
{
  data_allh <- rbind(data_allh, data_split[[SAMPLEID[i]]])
}

#===============================================================================
# Visualisation
#===============================================================================

plotUMAP <- function(x, title)
{
  ggplot() +
  geom_point(data = data, #reference data
             aes(x    = UMAP_1, 
                 y    = UMAP_2,
                 color = "A",
                 )
             ) +
  scale_color_manual(values = c("A" = "grey")) +
  guides(color = "none") +
  new_scale_color() +
  geom_point(data = x,
             aes(x     = UMAP_1, 
                 y     = UMAP_2,
                 color = celltype
                 )
             ) +
  labs(title = element_text(title)) +
  theme_void()
  }

#Embryonic cell type
p1 <- plotUMAP(data_ref, "Embrionic")

#Gastruloide cell type
p2 <- plotUMAP(data_allh, "Gastruloid") 

pdf(OUTPUTDIR1, height = 15, width = 30)
p1 + p2
dev.off()

#===============================================================================
# Visualisation of all cell types
#===============================================================================

all_celltypes <- unique(data_allh$celltype)
n             <- length(all_celltypes)

plot.scType.UMAP <- function(cell_type)
{
  ggplot() +
    geom_point(data = data, #all data
               aes(x    = UMAP_1, 
                   y    = UMAP_2,
                   color = "A",
               )
    ) +
    scale_color_manual(values = c("A" = "grey")) +
    guides(color = "none") +
    new_scale_color() +
    geom_point(data = data_allh,
               aes(x     = UMAP_1, 
                   y     = UMAP_2,
                   color = celltype == cell_type,
               )
    ) +
    scale_color_manual(values = c(alpha("grey",0), "red")) +
    guides(color = "none") +
    labs(title = element_text(cell_type)) +
    theme_void() +
    theme(title = element_text(size = 30))
}

plots_list <- list()
for(i in 1:n)
{
  p               <- plot.scType.UMAP(all_celltypes[i])
  plots_list[[i]] <- p
}

pdf(OUTPUTDIR2, height = 15*(ceiling(n/2)), width = 30)
plot_grid(plotlist = plots_list, ncol = 2)
dev.off()

#===============================================================================
# Visualisation of sample
#===============================================================================

all_samples <- unique(data_allh$sample)
n           <- length(all_samples)

plot.sample.UMAP <- function(cell_type)
{
  ggplot() +
    geom_point(data = data, #all data
               aes(x    = UMAP_1, 
                   y    = UMAP_2,
                   color = "A",
               )
    ) +
    scale_color_manual(values = c("A" = "grey")) +
    guides(color = "none") +
    new_scale_color() +
    geom_point(data = data_ref,
               aes(x     = UMAP_1, 
                   y     = UMAP_2,
                   color = sample == sample,
               )
    ) +
    scale_color_manual(values = c(alpha("grey",0), "red")) +
    guides(color = "none") +
    labs(title = element_text(cell_type)) +
    theme_void() +
    theme(title = element_text(size = 30))
}

plots_list <- list()
for(i in 1:n)
{
  p               <- plot.sample.UMAP(all_celltypes[i])
  plots_list[[i]] <- p
}

pdf(OUTPUTDIR3, height = 15*(ceiling(n/2)), width = 30)
plot_grid(plotlist = plots_list, ncol = 2)
dev.off()