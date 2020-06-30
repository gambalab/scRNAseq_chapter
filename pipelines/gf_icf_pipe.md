# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx.rds](https://drive.google.com/file/d/1TLYPcawqbtqApGDTXUZuwE9BJDpeUbPX/view?usp=sharing)
2. [tabulamuris_cell_info.rds](https://drive.google.com/file/d/1NjKix8SDCskJw1_IJtGib1IlSDdjfwOY/view?usp=sharing)

```R
# Install the package if necessary
if(!require(gficf))
{
  # Install required bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  BiocManager::install(setdiff(c("edgeR", "BiocParallel", "fgsea", "biomaRt","slingshot","tradeSeq"),rownames(installed.packages())),update = F)
  
  # install gficf package
  install.packages(pkgs = "gficf",repos = c("https://dibbelab.github.io/Rrepo/","https://cloud.r-project.org"))
}

# Load Tabula Muris dataset (already filterd)
M = readRDS("~/path/to/emptyDroplets_doublet_filtered_tabulamuris_mtx.rds")
cell.info = readRDS("~/path/to/tabulamuris_cell_info.rds")

# Gene filtering and data normalization
data = gficf::gficf(M = M,cell_proportion_min = .05,storeRaw = T,normalize = T,verbose = T)
rm(M);gc()

# PCA, data summarization
data = gficf::runPCA(data = data,dim = 50)

# Dimensionality reduction with UMAP 
data = gficf::runReduction(data = data,reduction = "umap",nt = 6,verbose = T,metric="manhattan")

# Identify clusters
data = gficf::clustcells(data = data,verbose = T,k = 50,nt = 6,community.algo = "louvian 2",resolution = .5,n.start = 50,n.iter = 250,dist.method = "manhattan")


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