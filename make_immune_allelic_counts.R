library(fishpond)
library(AnnotationHub)
library(ensembldb)
library(SummarizedExperiment)

ah <- AnnotationHub()
#query(ah, c("EnsDb","102","Mus musculus"))
edb <- ah[["AH89211"]]
txps <- transcripts(edb, return.type="DataFrame")
tx2gene <- txps[,c("tx_id","gene_id")]

cells <- paste0("cell",1:384)
metadata <- read.csv("metadata.csv")
stopifnot(all(metadata$cell == 1:384))
coldata <- data.frame(
  files = file.path("quants",cells,"quant.sf"),
  names = cells,
  cell = 1:384,
  cluster=metadata$cluster
)

cl <- c("T.Cells","Macrophages")
coldata <- coldata[coldata$cluster %in% cl,]
coldata <- coldata[order(coldata$cluster, coldata$cell),]
table(coldata$cluster)

se <- importAllelicCounts(coldata,
                          a1="alt", a2="ref",
                          format="wide",
                          tx2gene=tx2gene)

keep <- rowSums(assay(se) >= 10) >= 20
table(keep)
se <- se[keep,]

# add ranges
g <- genes(edb)
g <- g[rownames(se)]
rowRanges(se) <- g

save(se, file="allelic_counts_gene_filt.rda")

library(plyranges)
txps <- makeTx2Tss(edb, maxgap=50) %>%
  select(tx_id, gene_id, group_id, tss)

se <- importAllelicCounts(coldata,
                          a1="alt", a2="ref",
                          format="wide",
                          tx2gene=txps)

keep <- rowSums(assay(se) >= 10) >= 10
table(keep)
se <- se[keep,]

save(se, file="allelic_counts_tss_filt.rda")
