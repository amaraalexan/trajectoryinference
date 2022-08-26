# SlingShot 
# Amara Alexander 

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(grDevices)
library(GiNA)

setwd("")
# Data --------------
seurat_clustered <- readRDS("./seurat_clustered.rds")

# SlingShot -----
#convert to single cell experiment for input input into SlingShot

sce <- as.SingleCellExperiment(seurat_clustered, assay = "RNA")

sce <- slingshot(sce, reducedDim = 'PCA', clusterLabels = seurat_clustered@active.ident)
summary(sce$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

par(mar = c(1, 1, 1, 1))

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')


plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

seurat_clustered$pseudotime1 <- sce$slingPseudotime_1
seurat_clustered$pseudotime2 <- sce$slingPseudotime_2
FeaturePlot(seurat_clustered, c("pseudotime1", "pseudotime2"))

