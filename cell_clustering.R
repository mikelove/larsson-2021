suppressPackageStartupMessages(library(SummarizedExperiment))
load("total_counts_filt.rda")

assays(se) <- assays(se)[1] # just counts
metadata <- read.csv("metadata.csv")
stopifnot(all(se$cell == metadata$cell))
se$cluster <- metadata$cluster # add cluster from publication

library(scran)
library(scater)
sce <- as(se, "SingleCellExperiment")

# log scaled counts
set.seed(1)
clust <- quickCluster(sce, min.size=48)
sce <- computeSumFactors(sce, cluster=clust)
sce <- logNormCounts(sce)

# highly variable genes -> PCA
dec <- modelGeneVar(sce)
top <- getTopHVGs(dec, n=1000)
sce <- fixedPCA(sce, subset.row=top)
plotReducedDim(sce, dimred="PCA", color_by="cluster")

# TSNE similar to publication
sce <- runTSNE(sce, dimred="PCA")
plotReducedDim(sce, dimred="TSNE", color_by="cluster")

