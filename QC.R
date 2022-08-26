# SlingShot Workflow # 
#        7/17/22     #
#   Amara Alexander  #


# Set Up --------------
req.packages <- c("R.utils",
                  "scales",
                  "Seurat",
                  "stringr",
                  "tidyverse")


installed_packages <- req.packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(req.packages[!installed_packages], dependencies = TRUE, INSTALL_opts = "--no-lock", repos = "http://cran.us.r-project.org")
}
invisible(lapply(req.packages, library, character.only = TRUE))

options(stringsAsFactors = F)

#gunzip("./SRR12615659_wasp_Solo.out/SRR12615659_wasp_Solo.out/GeneFull/Features.stats.gz")
setwd()

## Reading in Data ------------------
#GeneFull/Filtered Data 
folder <- list.files("./pc_SoloOuts/data")

filelist <- lapply(folder, function(folder){
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})

pc_data <- merge(filelist[[1]], filelist[[2]], 
                       add.cell.ids = folder,
                       project= "prostate-scrna")

# Quality metrics -----------------
#determine mitochondrial ratio
pc_data$percent.mt <- PercentageFeatureSet(pc_data, pattern = "^MT-")
pc_data$percent.mt <- pc_data@meta.data$percent.mt/ 100

#Number of genes/UMI
pc_data$log10GenesPerUMI <- log10(pc_data$nFeature_RNA) / log10(pc_data$nCount_RNA)

# Create metadata dataframe
metadata <- pc_data@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)# Add cell IDs to metadata

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "[1]"))] <- "SRR12615659"
metadata$sample[which(str_detect(metadata$cells, "[2]"))] <- "SRR12615660"

#saving metadata 
pc_data@meta.data <- metadata

## Checking Quality Metrics -----
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
#Note: Sample 660 has 5000 more cells than sample 659

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
#Note: Cells in both samples have greater than 1000 UMIs (good: in-depth sequencing). Sample 660 has higher cell density.

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
#Note: Small shoulder could indicate biologically different types of cells

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 1500) +
  facet_wrap(~sample)
#Note: Some poor quality cells that can be filtered out

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
#Note: possibly more dying cells in sample 660, due to higher mitochondrial counts


#Pre-processing -----------------

##Filtering ------
filtered_seurat <- subset(x = pc_data, 
                          subset= (nUMI >= 1000) & 
                            (nGene >= 1500) & 
                            (log10GenesPerUMI > 0.75) & 
                            (percent.mt < 0.50))

# Filtering out genes that are expressed in less than 10 cells
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

### Confirming Metrics Used -----
# metadata_clean %>% 
#   ggplot(aes(x=sample, fill=sample)) + 
#   geom_bar() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells")
# 
# metadata_clean %>% 
#   ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   scale_x_log10() + 
#   theme_classic() +
#   ylab("Cell density") +
#   geom_vline(xintercept = 500)
# 
# metadata_clean %>% 
#   ggplot(aes(color=sample, x=nGene, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   theme_classic() +
#   scale_x_log10() + 
#   geom_vline(xintercept = 300)
# 
# metadata_clean %>% 
#   ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
#   geom_boxplot() + 
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells vs NGenes")
# 
# metadata_clean %>% 
#   ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
#   geom_point() + 
#   scale_colour_gradient(low = "gray90", high = "black") +
#   stat_smooth(method=lm) +
#   scale_x_log10() + 
#   scale_y_log10() + 
#   theme_classic() +
#   geom_vline(xintercept = 2000) +
#   geom_hline(yintercept = 1500) +
#   facet_wrap(~sample)
# 
# metadata_clean %>% 
#   ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
#   geom_density(alpha = 0.2) + 
#   scale_x_log10() + 
#   theme_classic() +
#   geom_vline(xintercept = 0.2)
#Note: Outliers have been removed with conservation of possibly dying cells with high mitochondrial percentage

# Create .RData object to load at any time
save(filtered_seurat, file="seurat_filtered.RData")

## Normalization  ----
#split seurat object by sample
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("SRR12615659", "SRR12615660")]

#peforming scTransform on both samples (normalizes, adjusts variance, and regresses out insignificant variation)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
}

# Save the split seurat object
saveRDS(split_seurat, "split_seurat.rds")

## Integration ----
#Tutorial Used: https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies 
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

PCAPlot(seurat_integrated,
        split.by = "sample")  

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)

# Save integrated seurat object
saveRDS(seurat_integrated, "integrated_seurat.rds")
