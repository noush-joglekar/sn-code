#/bin/R
# By Simon Hardwick, 2021
# Figure 2 for SnISOr-Seq paper

# Setup --------------------

library("dplyr")
library("ggplot2")
library("Seurat")

# Figure 2a --------------------

Gan62.data <- Read10X(data.dir = 'data/filtered_feature_bc_matrix')
Gan62 <- CreateSeuratObject(counts = Gan62.data, min.cells = 3, min.features = 200)
Gan62[["percent.mt"]] <- PercentageFeatureSet(Gan62, pattern = "^MT-")
head(Gan62@meta.data, 5)
VlnPlot(Gan62, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Gan62, features = "nCount_RNA")
VlnPlot(Gan62, features = "nFeature_RNA")

FeatureScatter(Gan62, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Gan62, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Gan62 <- subset(Gan62, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 4)
Gan62 <- NormalizeData(Gan62, normalization.method = "LogNormalize", scale.factor = 10000)
Gan62 <- FindVariableFeatures(Gan62, selection.method = "vst", nfeatures = 2000)
Gan62.top10 <- head(VariableFeatures(Gan62), 10)
Gan62.plot1 <- VariableFeaturePlot(Gan62)
Gan62.plot2 <- LabelPoints(plot = Gan62.plot1, points = Gan62.top10, repel = TRUE)
Gan62.plot2
Gan62.all.genes <- rownames(Gan62)
Gan62 <- ScaleData(Gan62, features = Gan62.all.genes)

Gan62 <- RunPCA(Gan62, features = VariableFeatures(object = Gan62))
print(Gan62[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Gan62, dims = 1:2, reduction = "pca")
DimPlot(Gan62, reduction = "pca")
DimHeatmap(Gan62, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Gan62, dims = 1:15, cells = 500, balanced = TRUE)
Gan62 <- JackStraw(Gan62, num.replicate = 100)
Gan62 <- ScoreJackStraw(Gan62, dims = 1:20)
JackStrawPlot(Gan62, dims = 1:15)
ElbowPlot(Gan62)
Gan62 <- FindNeighbors(Gan62, dims = 1:13)
Gan62 <- FindClusters(Gan62, resolution = 0.6)
Gan62 <- RunUMAP(Gan62, dims = 1:13)
DimPlot(Gan62, reduction = "umap", label = TRUE)

# Figure 2b generated with Prism using file 'data/snisor_readMetrics.csv'

# Figure 2c --------------------

gene_corr <- read.csv('data/Gan62vs64_correlation_tpms.csv', header = TRUE)

ggplot(data = gene_corr, aes(x=log(Gan62), y=log(Gan64))) + geom_point()

ggplot(data = gene_corr, aes(x=log10(Gan62+1), y=log10(Gan64+1))) + geom_point(shape=1, size=2) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cor.test(log(gene_corr[,2]+1),log(gene_corr[,3]+1))

# Figure 2d --------------------

gene_exp <- read.csv('data/Gan62_combined_tpms.csv', header = TRUE)

ggplot(data = gene_exp, aes(x=log10(control+1), y=log10(exome+1))) + geom_point(shape=1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cor.test(log(gene_exp[,2]+1),log(gene_exp[,3]+1))
