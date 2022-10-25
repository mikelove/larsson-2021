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
coldata <- data.frame(
  files = file.path("quants",cells,"quant.sf"),
  names = cells,
  cell = 1:384
)
se <- importAllelicCounts(coldata,
                          a1="alt", a2="ref",
                          format="assays",
                          tx2gene=tx2gene,
                          dropInfReps=TRUE)

assay(se, "counts") <- assay(se, "a1-counts") + assay(se, "a2-counts")
assays(se) <- assays(se)[c(7, 1:6)]
keep <- rowSums(assay(se) >= 10) >= 48
table(keep)
se <- se[keep,]

# add ranges
g <- genes(edb)
g <- g[rownames(se)]
rowRanges(se) <- g

save(se, file="total_counts_filt.rda")
