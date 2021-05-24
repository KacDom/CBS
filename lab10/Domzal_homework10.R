library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


pbmc.data <- Read10X(data.dir = "./filtered_matrices_mex/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
pbmc <- ScaleData(pbmc)

"1. Analyze this dataset in the same way as Seurat's guide with PBMC3k. 
Apply PCA, clustering, and t-SNE to reproduce Figure 3b on Zheng et al. 2017."
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
pbmc <- FindNeighbors(pbmc, graph.name = "test", dims = 1:12)
pbmc <- FindClusters(pbmc, graph.name = "test", resolution = 0.148)
pbmc <- RunTSNE(pbmc, dims = 1:12)
DimPlot(pbmc, reduction = "tsne", label = TRUE, label.size = 7) + ggtitle("tSNE of cells clutered by FindClusters")



"2. Create a hierachical clustering by applying K-means clustering (or FindClusters) to cells defined
by each of 10 cluster. Try to find a suitable number of clusters (k) for each sub-population.
Present overall t-SNE visualization with all clusters - make sure that hierarchy (group) is visualized.
visualize t-SNE for each of 10 clusters, with their sub-clusters."
num_of_clusters = 10
pbmc$sub_cluster <- as.character(Idents(pbmc))
plot_list = list()

pdf("Problem2_sub.pdf")
for (idx in 0:(num_of_clusters-1)){
  cluster <- subset(pbmc, idents = idx)
  cluster <- FindNeighbors(cluster, dims = 1:12)
  cluster <- FindClusters(cluster, resolution = 0.2)
  pbmc$sub_cluster[Cells(cluster)] <- paste(idx, Idents(cluster), sep = '_')
  print(DimPlot(subset(pbmc, idents = idx), group.by = "sub_cluster", label = TRUE, label.size = 5) + ggtitle(paste("tSNE of subcluster number: ", idx)))
}
dev.off()

DimPlot(pbmc, group.by = "sub_cluster", label = TRUE, label.size = 5) + ggtitle("tSNE of cells with subclusters")


