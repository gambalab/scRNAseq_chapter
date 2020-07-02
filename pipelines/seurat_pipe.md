# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx.rds](https://drive.google.com/file/d/1TLYPcawqbtqApGDTXUZuwE9BJDpeUbPX/view?usp=sharing)
2. [tabulamuris_cell_info.rds](https://drive.google.com/file/d/1NjKix8SDCskJw1_IJtGib1IlSDdjfwOY/view?usp=sharing)


# Tabule Muris Example

## 1. First we install all the necessary packages.
```R
# Install the package if necessary
if(!require(Seurat))
{
  install.packages(pkgs = "Seurat")
  library(Seurat)
}
```
## 2. We load into the R environment the downloaded Tabula Muris datasets

```R
# Load Tabula Muris dataset (already filterd)
tabData=readRDS('/path/to/emptyDroplets_doublet_filtered_tabulamuris_mtx.rds')
info=readRDS("/path/to/tabulamuris_cell_info.rds")

# We also filter genes in less than 5% of total cells
tabmur <- CreateSeuratObject(counts = tabData, project = "tabulamuris", min.cells = round(dim(tabData)[2]*5/100), min.features = 0, meta.data=info)
```

## 3. Data normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# uses log transformation of the CPM to reduce cell depth variability.
tabmur <- NormalizeData(tabmur, normalization.method = "LogNormalize", scale.factor = 10000)
rm(tabData);gc()
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties. Dimensionality reduction methods are essential for clustering, visualization, and summarization of scRNA-seq data. PCA is used to summarise a dataset throughout the top N principal components. The number of PCA to use is usually determined by manually inspecting the elbow plot (*ElbowPlot* function), in which principal components are plotted as a function of the variability they account for, and the number of PCA to use is determined by the point in which an elbow is observed. Additional methods can be used, including jackstraw 
```R
# Feature selection with top 2,000 High Variable Genes (HGV)
tabmur <- FindVariableFeatures(object = tabmur, selection.method="vst", nfeatures = 2000)
top10 <- head(VariableFeatures(tabmur), 10)

# plot HGV
plot1 <- VariableFeaturePlot(tabmur)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Seurat suggest to rescale expression before PCA.
# Shifts the expression of each gene, so that the mean expression across cells is 0 setting do.center=TRUE
# Scales the expression of each gene, so that the variance across cells is 1 setting do.scale=TRUE
tabmur <- ScaleData(tabmur, do.center = TRUE, do.scale = TRUE)

# Finally we run PCA
tabmur <- RunPCA(tabmur, features = VariableFeatures(object = tabmur),npcs = 50)
ElbowPlot(tabmur,ndims=50)

# Dimensionality reduction with UMAP 
tabmur <- RunUMAP(tabmur, dims = 1:50)
```


## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity.
**Seurat** constructs a KNN graph based on the euclidean distance in PCA space, and refines the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the *FindNeighbors* function, and takes as input the first 50 PCs.

To cluster the cells, **Seurat** applys a modularity optimization technique such as the Louvain algorithm (default) or SLM. The *FindClusters* function implements this procedure. When running a graph-based clustering, it is necessary to set the resolution parameter for the community detection algorithm based on modularity optimization. The resolution parameter is correlated to the scale of observing communities. In particular, the higher is the resolution parameter, the larger is the number of smaller communities. 
Here, we set the resolution parameter to 0.5

```R
# Identify clusters
tabmur <- FindNeighbors(tabmur, dims = 1:50)
tabmur <- FindClusters(tabmur, resolution = 0.5)
```

## 6. Plots
```R
DimPlot(tabmur, reduction = "umap")
```


