# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx.rds](https://drive.google.com/file/d/1TLYPcawqbtqApGDTXUZuwE9BJDpeUbPX/view?usp=sharing)
2. [tabulamuris_cell_info.rds](https://drive.google.com/file/d/1NjKix8SDCskJw1_IJtGib1IlSDdjfwOY/view?usp=sharing)


# Tabule Muris Example

## 1. First we install all the necessary packages.
```R
# Install the package if necessary
if(!require(gficf))
{
  # Install required bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  BiocManager::install(setdiff(c("edgeR", "BiocParallel", "fgsea", "biomaRt","slingshot","tradeSeq"),rownames(installed.packages())),update = F)
  
  # install gficf package
  install.packages(pkgs = "gficf",repos = c("https://dibbelab.github.io/Rrepo/","https://cloud.r-project.org"))
  
  library(gficf)
}
```

## 2. We load into the R environment the downloaded Tabula Muris datasets

```R
# Load Tabula Muris dataset (already filterd)
M = readRDS("/path/to/emptyDroplets_doublet_filtered_tabulamuris_mtx.rds")
cell.info = readRDS("/path/to/tabulamuris_cell_info.rds")
```

## 3. Gene Filtering and normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# Gene filtering and data normalization
data = gficf::gficf(M = M,cell_proportion_min = .05,storeRaw = T,normalize = T,verbose = T)
rm(M);gc()
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```R
# PCA for data summarization
data = gficf::runPCA(data = data,dim = 50)

# Dimensionality reduction with UMAP 
data = gficf::runReduction(data = data,reduction = "umap",nt = 6,verbose = T,metric="manhattan")
```

## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.

```R
# Identify clusters
data = gficf::clustcells(data = data,verbose = T,k = 50,nt = 6,community.algo = "louvian 2",resolution = .5,n.start = 50,n.iter = 250,dist.method = "manhattan")
```

## 6. Plots
```R
# add info to cells
data$embedded$cell.id = rownames(data$embedded)
data$embedded$tissue = as.character(cell.info$tissue[match(data$embedded$cell.id,cell.info$id)])
data$embedded$cell_ontology_class = as.character(cell.info$cell_ontology_class[match(data$embedded$cell.id,cell.info$id)])

# Plot cells by ontology
p1 = gficf::plotCells(data,colorBy = "cell_ontology_class",pointSize = .1) + xlab("UMAP 1") + ylab("UMAP 2")
print(p1)

# Plot Cells by Clusters
p2 = gficf::plotCells(data,colorBy = "cluster",pointSize = .1) + xlab("UMAP 1") + ylab("UMAP 2")
print(p2)
```

|        Plot p1 (by ontology)     |       Plot p2 (by clusters)     |
|-----------------------------------|---------------------------------|
|![GFICF_umap_ontology.png](https://github.com/gambalab/scRNAseq_chapter/blob/master/plots/GFICF_umap_ontology.png?raw=true)|![GFICF_umap_clusters.png](https://github.com/gambalab/scRNAseq_chapter/blob/master/plots/GFICF_umap_clusters.png?raw=true)|



