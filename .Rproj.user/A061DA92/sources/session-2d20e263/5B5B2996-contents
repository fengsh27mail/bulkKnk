# HGC: fast hierarchical clustering for large-scale single-cell data
## Introduction

`HGC` (short for Hierarchical Graph-based Clustering) is an R package for 
conducting hierarchical clustering on large-scale single-cell RNA-seq 
(scRNA-seq) data. The key idea is to construct a dendrogram of cells on 
their shared nearest neighbor (SNN) graph. `HGC` provides functions for 
building cell graphs and for conducting hierarchical clustering on the graph. 
Experiments on benchmark datasets showed that `HGC` can reveal the 
hierarchical structure underlying the data, achieve state-of-the-art 
clustering accuracy and has better scalability to large single-cell data. 
For more information, please refer to the paper on 
[bioinformatics](https://doi.org/10.1093/bioinformatics/btab420) 
or the preprint of `HGC` on 
[bioRxiv](https://doi.org/10.1101/2021.02.07.430106).

## Installation

`HGC` has been published on 
[bioconductor](https://bioconductor.org/packages/HGC/).

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HGC")
```

`HGC` could also be installed from Github.

```{r Github install, eval = FALSE}
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("XuegongLab/HGC")
```

Different branches here provide variants of `HGC` for convenience. The
`HGC` packages from `bioconductor` and Github master branch are built in 
R 4.1. For the users with lower R versions, we suggest to use the `HGC` in 
[HGC4oldRVersion](https://github.com/XuegongLab/HGC/tree/HGC4oldRVersion) 
branch.

```{r Github install, eval = FALSE}
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("XuegongLab/HGC", ref = "HGC4oldRVersion")
```

For the users just interested in the core hierarchical clustering functions, 
they could reference the 
[HGC_core](https://github.com/XuegongLab/HGC/tree/HGC_core) branch.

```{r Github install, eval = FALSE}
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("XuegongLab/HGC", ref = "HGC_core")
```

## Quick Start

### Input data
`HGC` takes a matrix as input where row represents cells and column 
represents features. Preprocessing steps like normalization and dimension 
reduction are necessary so that the constructed graph can capture the 
manifold underlying the single-cell data. We recommend users to follow 
the standard preprocessing steps in 
[`Seurat`](https://satijalab.org/seurat/articles/get_started.html). 
As a demo input, we stored the top 25 principal components of the 
Pollen dataset ([Pollen et al.](https://www.nature.com/articles/nbt.2967)) 
in `HGC`. The dataset contains 301 cells with two known labels: labels at 
the tissue level and the cell line level.

```{r, message=FALSE, warning=FALSE}
library(HGC)

data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
Pollen.Label.Tissue <- Pollen[["Tissue"]]
Pollen.Label.CellLine <- Pollen[["CellLine"]]

dim(Pollen.PCs)
table(Pollen.Label.Tissue)
table(Pollen.Label.CellLine)
```

### Run HGC

There are two major steps for conducting the hierarchical clustering 
with `HGC`: the graph construction step and the dendrogram construction 
step. `HGC` provides functions for
building a group of graphs, including the k-nearest neighbor graph (KNN), 
the shared nearest neighbor graph (SNN), the continuous k-nearest neighbor 
graph (CKNN), etc. These graphs are saved as `dgCMatrix` supported by 
R package `Matrix`. Then `HGC` can directly build a hierarchical tree 
on the graph. A self-built graph or graphs from other pipelines stored 
as `dgCMatrix` are also supported.

```{r}
Pollen.SNN <- SNN.Construction(mat = Pollen.PCs, k = 25, threshold = 0.15)
Pollen.ClusteringTree <- HGC.dendrogram(G = Pollen.SNN)
```

The user could also give `HGC.dendrogram` an adjacency matrix directly, please 
reference to check the accepted data structures in the function documentation.  
For instance, read a matrix from `igraph` object and use it to run `HGC`.

```{r}
require(igraph)
g <- sample_gnp(10, 2/10)
G.ClusteringTree <- HGC.dendrogram(G = g)
```

The output of `HGC` is a standard tree following the data structure `hclust()` 
in R package `stats`. The tree can be cut into specific number of clusters 
with the function `cutree`.  

```{r}
cluster.k5 <- cutree(Pollen.ClusteringTree, k = 5)
```

### Visualization

With various published methods in R, results of `HGC` can be visualized easily. 
Here we use the R package `dendextend` as an example to visualize the results 
on the Pollen dataset. The tree has been cut into five clusters. And for a 
better visualization, the height of the tree has been log-transformed.

```{r, fig.height = 4.5}
Pollen.ClusteringTree$height = log(Pollen.ClusteringTree$height + 1)
Pollen.ClusteringTree$height = log(Pollen.ClusteringTree$height + 1)

HGC.PlotDendrogram(tree = Pollen.ClusteringTree,
                    k = 5, plot.label = FALSE)
```
We can also add a colour bar of the known label under the dendrogram as a 
comparison of the achieved clustering results.

```{r, fig.height = 4.5}
Pollen.labels <- data.frame(Tissue = Pollen.Label.Tissue,
                            CellLine = Pollen.Label.CellLine)
HGC.PlotDendrogram(tree = Pollen.ClusteringTree,
                    k = 5, plot.label = TRUE, 
                    labels = Pollen.labels)
```

### Evaluation of the clustering results

For datasets with known labels, the clustering results of `HGC` can be 
evaluated by comparing the consistence between the known labels and the 
achieved clusters. Adjusted Rand Index (ARI) is a wildly used statistics 
for this purpose. Here we calculate the ARIs of the clustering results at 
different levels of the dendrogram with the two known labels. 

```{r}
ARI.mat <- HGC.PlotARIs(tree = Pollen.ClusteringTree,
                        labels = Pollen.labels)
```

### Heatmap with the HGC tree

With the help of `pheatmap` package, we can combine the 
`HGC` clustering tree with the heatmap of gene expression 
data or low-dimensional data.

```{r}
# Input the clustering tree to pheatmap function
require(pheatmap)
pheatmap(mat = Pollen.PCs, cluster_rows = Pollen.ClusteringTree, 
        cluster_cols = FALSE, show_rownames = FALSE)
```

## Time complexity analysis of HGC

Our work shows that the dendrogram construction in `HGC` has a linear time 
complexity. For advanced users, `HGC` provides functions to conduct time 
complexity analysis on their own data. The construction of the dendrogram 
is a recursive procedure of two steps: 1. find the nearest neighbour pair, 
2. merge the node pair and update the graph. For different data structures of 
graph, there's a trade-off between the time consumptions of the two steps. 
Generally speaking, storing more information about the graph makes it faster 
to find the nearest neighbour pair (step 1) but slower to update the graph 
(step 2). We have experimented several datasets and chosen the best data 
structure for the overall efficiency. 

The key parameters related to the time consumptions of the two steps are the 
length of the nearest neighbor chains and the number of nodes needed to be 
updated in each iteration, respectively (for more details, please refer to 
our [preprint](https://doi.org/10.1101/2021.02.07.430106)).`HGC` provides 
functions to record and visualize these parameters.

```{r}
Pollen.ParameterRecord <- HGC.parameter(G = Pollen.SNN)

HGC.PlotParameter(Pollen.ParameterRecord, parameter = "CL")
HGC.PlotParameter(Pollen.ParameterRecord, parameter = "ANN")
```
