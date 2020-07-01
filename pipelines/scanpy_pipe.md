# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx_transposed.mtx](https://drive.google.com/file/d/16G9Gcojd6CqoK_VNSiz2i0vMFCmd9Q9V/view?usp=sharing)
2. [Tabulamuris_cmd.txt](https://drive.google.com/file/d/1ngJ45fOzgY6pm9gNZdCKjXdG7lmJmY03/view?usp=sharing)
3. [Tabulamuris_genes.txt](https://drive.google.com/file/d/1L3IGc59iVLwT2HE7sGBu3R_DhzjxT7oF/view?usp=sharing)


# Tabule Muris Example

## 1. First we install all the necessary packages.
Please, follow the istruction present on the official page [HERE](https://scanpy.readthedocs.io/en/stable/installation.html)


## 2. We load into the Python environment the downloaded Tabula Muris datasets

```python
(python environment)

import numpy as np
import pandas as pd
import anndata as ad
print(ad.__version__)
import scanpy as sc
print(sc.__version__)

inp_dir = "/path/to/directory/with/data/"

#file to save the final object in a python-compatible format
results_file = ''.join([inp_dir, "ExprMatrix.h5ad"])

#data frame containing information about cells
obs=pd.read_csv(''.join([inp_dir, "Tabulamuris_cmd.txt"]),sep="\t", index_col=0) 

#data frame containing genes annotation
var=pd.read_csv(''.join([inp_dir, "Tabulamuris_genes.txt"]), index_col=0, sep="\t")

#get the transposed matrix filtered for empty droplets and doublets
adata=ad.read_mtx(''.join([inp_dir, "emptyDroplets_doublet_filtered_tabulamuris_mtx_transposed.mtx"]))

#add metadata to the scanpy object
adata.obs=obs
adata.var=var
```


## 3. Gene Filtering and normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```python
#Scanpy suggests to filter out genes wich are expressed in less than 3 cells.

sc.pp.filter_genes(adata, min_cells=3)

#Scanpy provides a normalization that generates CPM values and it performs logarithm on pseudocounts.
#This data is stored in the raw data compartment of adata.
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```python
# Gene filtering based on Hiigh Variable Genes (HVGs)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

#plot the metrics used for the selection of HVG
sc.pl.highly_variable_genes(adata, save='.pdf') 

#check the number of highly variable genes
sum(adata.var['highly_variable']) 

#Before PCA reduction, Scanpy allows to regress out the influence that certain variables could have on data before reducing dimensions.
#It is suggested to use it on the percentage of mitochondrial gene expression and on the total number of reads per cell (not used for the sake of this review).
#Then data can be scaled (not applied) and the dimension is reduced to the PCA components (first 50 by default).
sc.tl.pca(adata, svd_solver='arpack')

#elbow plot to observe the relationship between PCs and variance 
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save=".pdf")

# UMAP embedding
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)
```

## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.

```python
sc.tl.louvain(adata, resolution=0.5)
adata.write(results_file)
```