# Clustering 

library(Seurat)
library(dplyr)
# Data --------------
seurat_integrated <- readRDS("./integrated_seurat.rds")

# Cluster------

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Explore resolutions
seurat_integrated@meta.data %>% 
  View()
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
#14 clusters

get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(c(0:10), get_conserved)

#Find conserved markers for cluster 11, 12, & 13 separately due to size 
cluster_11 <- FindConservedMarkers(seurat_integrated, 
                                      ident.1 = 11,
                                      grouping.var = 'sample',
                                      only.pos = TRUE)
cluster_12 <- FindConservedMarkers(seurat_integrated, 
                                   ident.1 = 12,
                                   grouping.var = 'sample',
                                   only.pos = TRUE)
cluster_13 <- FindConservedMarkers(seurat_integrated, 
                                   ident.1 = 13,
                                   grouping.var = 'sample',
                                   only.pos = TRUE)
# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (SRR12615659_avg_log2FC + SRR12615660_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

#Assigning Clusters
# Human Protein Atlas used for assigning clusters (proteinatlas.org) 
seurat_integrated<- RenameIdents(seurat_integrated, "0"  = "Prostatic Glandular", "1" = "Smooth Muscle", "2" = "Endothelial", "3" = "Endothelial","4"= "Basal Prostatic", "5"= "Endothelial", "6" = "Urothelial", "7" = "Smooth Muscle", "8" ="T-cell", "9"="Macrophage","10"= "Urothelial","11" ="Fibroblast", "12"= "Prostatic Glandular", "13"= "T-cell")

DimPlot(seurat_integrated, reduction = "umap", label= T, pt.size = 0.3)

saveRDS(seurat_integrated, file = "seurat_clustered.rds")
