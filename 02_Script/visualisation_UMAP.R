#===============================================================================
#Library
#===============================================================================

library(ggplot2)
library(ggnewscale)
library(Seurat)

#===============================================================================
#Data
#===============================================================================

INPUTDIR   <- snakemake@input[["input"]]
OUTPUTDIR  <- snakemake@output[["output"]]
integrated <- readRDS(INPUTDIR)

SAMPLEID   <- c("72h", "80h", "86h", "96h", "mouse gastrulation")

#Create a data frame 4 x number of cells, with each cells assigned the 
#parameters UMAP_1, UMAP_2, origin, celltype
data          <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
data$origin   <- unlist(as.list(integrated@meta.data[["orig.ident"]]))
data$celltype <- unlist(as.list(integrated@meta.data[["celltype"]]))

#Split for each condition
data_split    <- split(data, data$origin)
n             <- length(data_split)

data_ref  <- data_split[[SAMPLEID[n]]]  
data_allh <- data_split[[SAMPLEID[1]]]  
for(i in 2:(n-1))
{
  data_allh <- rbind(data_allh, data_split[[SAMPLEID[i]]])
}

integrated <- NA

#===============================================================================
#Visualisation
#===============================================================================

plotUMAP <- function(x, title)
{
  ggplot() +
  geom_point(data = data,
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

pdf(OUTPUTDIR, height = 15, width = 30)
p1 + p2
dev.off()

