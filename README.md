# Single-cell RNA sequencing analysis: a step by step overview
#### S. Slovin<sup>1,+</sup>, A. Carissimo<sup>2,+</sup>, F. Panariello<sup>1,+</sup>, A. Grimaldi<sup>1</sup>, V. Bouché<sup>1</sup>, G. Gambardella<sup>1,3,++</sup>, D. Cacchiarelli<sup>1,4,++</sup>

<sup>1</sup>Telethon Institute of Genetics and Medicine (TIGEM), Armenise/Harvard Laboratory of Integrative Genomics, Pozzuoli, Italy.  
<sup>2</sup>Istituto per le Applicazioni del Calcolo "Mauro Picone," Consiglio Nazionale delle Ricerche, Naples, Italy.  
<sup>3</sup>Department of Chemical Materials and Industrial Engineering, University of Naples “Federico II”, Naples, Italy.  
<sup>4</sup>Department of Translational Medicine, University of Naples "Federico II", Naples, Italy.

## Abstract
Thanks to innovative sample-preparation and sequencing technologies, gene expression in individual cells can now be measured for thousands of cells in a single experiment. Since its introduction, single-cell RNA sequencing (scRNA-seq) approaches have revolutionized the genomics field as they created unprecedented opportunities for resolving cell heterogeneity by exploring gene expression profiles at a single-cell resolution. However, the rapidly evolving field of scRNA-seq invoked the emergence of various analytics approaches aimed to maximize the full potential of this novel strategy. Unlike population-based RNA sequencing approaches, scRNA seq necessitates comprehensive computational tools to address high data complexity and keep up with the emerging single-cell associated challenges. Despite the vast number of analytical methods, a universal standardization is lacking. While this reflects the fields' immaturity, it may also encumber a newcomer to blend in.

In this review, we aim to bridge over the above-mentioned hurdle and propose four ready-to-use pipelines for scRNA-seq analysis easily accessible by a newcomer, that could fit various biological data types. Here we provide an overview of the currently available single-cell technologies for cell isolation and library preparation and a step by step guide that covers the entire canonical analytic workflow to analyse scRNA-seq data including read mapping, quality controls, gene expression quantification, normalization, feature selection, dimensionality reduction, and cell clustering useful for trajectory inference and differential expression. Such workflow guidelines will escort novices as well as expert users in the analysis of complex scRNA-seq datasets, thus further expanding the research potential of single-cell approaches in basic science, and envisaging its future implementation as best practice in the field.

The full article [(Slovin et al. 2020)](https://) will be soon available ...

## Pipelines Examples
1. [Monocle 3](https://github.com/gambalab/scRNAseq_chapter/blob/master/pipelines/monocle3_pipe.md) (R) 
2. [GF-ICF](https://github.com/gambalab/scRNAseq_chapter/blob/master/pipelines/gf_icf_pipe.md) (R)
3. Seurat (R)
4. [Scanpy](https://github.com/gambalab/scRNAseq_chapter/blob/master/pipelines/scanpy_pipe.md) (Python)

## Original Pubblications of used tools:
Monocle [(Cacchiarelli et al. 2014)](https://www.nature.com/articles/nbt.2859) 
GF-ICF [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract) 
Seurat [(Stuart et al. 2019)](https://doi.org/10.1016/j.cell.2019.05.031) 
Scanpy [(Wolf et al. 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0) 

## Contents of the article

1. Introduction / Background

2. The laboratory workflow of scRNA-seq
3. The computational workflow of scRNA-seq
4. Raw reads Demultiplexing, Alignment and Expression Quantification 
    1. Demultiplexing
    2. Mapping and expression quantification
5. Quality Control and cell filtering: How to identify viable cells
    1. Identify empty barcoded drops
    2. Multiplets identification
    3. Cells lysis
6. Start working with the Scanpy, Seurat, Monocle and gf-icf's pipelines:
7. Gene filtering: how to remove “noisy” genes
8. Data normalization: how to make gene expression comparable across individual cells
9. Feature selection: how to discard “uninformative” genes
10. Dimensionality reduction: how to summarize and visualize scRNA-seq data
    1. Linear dimensional reduction: for the summarization of scRNA-seq data
    2. Nonlinear dimensionality reduction for the visualization of scRNA-seq data
11. Clustering Analysis: how to identify cellular sub-populations
12. Differential Expression: how to annotate cell populations
13. Results evaluation and comparison among the implemented pipelines
14. Additional analyses: how to reconstruct cell transcriptional dynamics
15. Discussion and future directions

## Figures

![alt text](https://github.com/gambalab/scRNAseq_chapter/blob/master/figures_HiRes/Figure_01.png?raw=true)

<b>Figure 1 - Single-cell RNA sequencing workflow.</b> The scRNA-seq procedure consists of six key steps. (i) Samples are dissociated into a single-cell suspension. (ii) As lysed cells might bias the data and cause high noise interference, it is essential to maximize the quality of the input material and assess cell viability. (iii) If the viability is lower than 90%, dead cells should be filtered either by centrifugation (i.e density gradient) or immunodepletion (i.e FACS or magnetic sorting). (iv) Single cells are captured and isolated in different ways, depending on the technique of choice. Microfluidics-based scRNA-seq technologies encapsulate single cells within water-in-oil droplets together with unique primers attached to microparticle surface and lysis buffer. Then, each lysed cell's mRNA content is captured by the poly-A tail domain of a single primer and labeled with UMI and cell-specific barcodes. Several errors can occur during this step, like multiple cells or microparticles captured in a single droplet (i.e. multiplets), and sub-Poisson loading trade-offs, such as empty barcoded drops . (v) Captured mRNA transcripts from droplets are then collected, reverse-transcribed, and (vi) amplified in pools to be used for standard sequencing platforms. During library construction, cDNA molecules are tagged with sample-specific indexes allowing multiplexing of different captures in the same sequencing run. Further computational demultiplexing will use such barcode information to sort samples, cells, and transcripts.
<hr/>

![alt text](https://github.com/gambalab/scRNAseq_chapter/blob/master/figures_HiRes/Figure_02.png?raw=true)
<b>Figure 2 - Computational analysis of single cell RNA sequencing.</b> ScRNA-seq analysis embraces six underlying steps, including raw-data pre-processing, filtering via QC covariants, normalization, feature selection, linear dimensional reduction, visualization, and clustering: (i) Raw reads are processed and quantified to generate gene/barcode matrices. (ii) Cells in the count matrix are then filtered to avoid misinterpretation of ambient gene expression, apoptotic cells, and multiplets. (iii) Count reads normalization is required, as the analysis is disrupted by low input and weak SNR, following which data is primed for downstream analysis. (iv) A lesser number of highly variable features are selected for the purpose of realizing a faster and accurate procedure. (v) Based on the designated genes, a PCA is performed to lower data dimensionality. (vi) Clustering and nonlinear dimensionality reduction steps utilize a subset of significant principal components to overcome data noisiness. Subsequently, cells are clustered and visualized based on their PCA scores.
<hr/>

![alt text](https://github.com/gambalab/scRNAseq_chapter/blob/master/figures_HiRes/Figure_03.png?raw=true)
<b>Figure 3 - Cell QC on Tabula Muris dataset.</b> (A) Detection of empty droplet by using emptyDrop function from DropUtils R package on Tabula Muris dataset. (B) Identification of cell multiplets in each independent run of Tabula Muris Dataset. (C) Distribution number of detected genes across the cells in the Tabula muris dataset. (D) PCA components as a function of their percentage of explained variance on Tabula Muris dataset (elbow plot).
<hr/>

![alt text](https://github.com/gambalab/scRNAseq_chapter/blob/master/figures_HiRes/Figure_04.png?raw=true)
<b>Figure 4 - Cell clustering, comparison and evaluation of the implemented pipelines.</b> (A-D) UMAP visualization produced by each one of the 4 implemented pipelines where cells are colored according to the cluster in which they fall. (E) Agreement of identified cell clusters among the four implemented pipelines. (F) Cluster results of each independent pipeline is used to hierarchically cluster cell types and reconstruct cell lineage.
