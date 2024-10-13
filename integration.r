library(dplyr)
library(patchwork)
library(Seurat)
library(httpgd)
library(hdf5r)

data_0_hpa_B <- Read10X(data.dir = "/Users/ladislavbucek/Desktop/UK/Mgr./labgenexp/ravi_data_tutorial/data/raw/regeneration/STAR/0_hpa_B/filtered")
data_1_hpa_B <- Read10X(data.dir = "/Users/ladislavbucek/Desktop/UK/Mgr./labgenexp/ravi_data_tutorial/data/raw/regeneration/STAR/1_hpa_B/filtered")
data_3_hpa_A <- Read10X(data.dir = "/Users/ladislavbucek/Desktop/UK/Mgr./labgenexp/ravi_data_tutorial/data/raw/regeneration/STAR/3_hpa_A/filtered")
data_12_hpa_A <- Read10X(data.dir = "/Users/ladislavbucek/Desktop/UK/Mgr./labgenexp/ravi_data_tutorial/data/raw/regeneration/STAR/12_hpa_A/filtered")

pbmc.data <- Read10X(data.dir = "/Users/ladislavbucek/Desktop/UK/Mgr./labgenexp/tutorial seurat/introductory vignettes/data/raw/filtered_gene_bc_matrices/hg19")


head(data_0_hpa_B)
head(pbmc.data)


obj_0_hpa_B <- CreateSeuratObject(counts = data_0_hpa_B, project = "ravi_0_hpa_B", min.cells = 3, min.features = 200)
obj_0_hpa_B[["percent.mt"]] <- PercentageFeatureSet(obj_0_hpa_B, pattern = "^mt-")


# Visualize QC metrics as a violin plot
VlnPlot(obj_0_hpa_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(obj_0_hpa_B, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj_0_hpa_B, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
