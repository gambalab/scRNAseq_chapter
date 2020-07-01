# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx.rds](https://drive.google.com/file/d/1TLYPcawqbtqApGDTXUZuwE9BJDpeUbPX/view?usp=sharing)
2. [Tabulamuris_cmd.txt](https://drive.google.com/file/d/1ngJ45fOzgY6pm9gNZdCKjXdG7lmJmY03/view?usp=sharing)
3. [Tabulamuris_genes.txt](https://drive.google.com/file/d/1L3IGc59iVLwT2HE7sGBu3R_DhzjxT7oF/view?usp=sharing)


# Tabule Muris Example

## 1. First we install all the necessary packages.
Additional information on possible errors during installation process can be found [HERE](https://cole-trapnell-lab.github.io/monocle3/docs/installation/)

```R
if(!require(monocle3))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(version = "3.10")
  
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
  
  install.packages("devtools")
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
  
  devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
  
  library(monocle3)
}
```


## 2. We load into the R environment the downloaded Tabula Muris datasets

```R
# Load Tabula Muris dataset (already filterd)
mtx=readRDS("/path/to/emptyDroplets_doublet_filtered_tabulamuris_mtx.rds") #get the matrix filtered for empty droplets and doublets
cell_metadata=read.table("/path/to/Tabulamuris_cmd.txt", header=T,row.names=1) #data frame containing information about cells
gene_metadata=read.table("/path/to/Tabulamuris_genes.txt", header=T, row.names=1) #data frame containing genes annotation
```

## 3. Gene Filtering and normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# Gene filtering and data normalization
cds <- new_cell_data_set(as(mtx, "sparseMatrix"),cell_metadata = cell_metadata,gene_metadata = gene_metadata)
rm(mtx);gc()
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```R
# PCA for data summarization
cds <- preprocess_cds(cds, num_dim = 50)

# Dimensionality reduction with UMAP 
cds <- reduce_dimension(cds) #dimensionality reduction, default value is UMAP
```


## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.

```R
# Identify clusters
cds = cluster_cells(cds, cluster_method="louvain") #cell clustering using louvain algorithm
```

## 6. Plots
```R
plot_cells(cds)
```


