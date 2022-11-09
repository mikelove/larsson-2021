suppressPackageStartupMessages(library(SummarizedExperiment))

load("allelic_counts_gene_filt.rda")
se$cluster <- factor(se$cluster)
levels(se$cluster) <- c("macrophage","T-cell")

se <- se[,order(se$cluster,se$cell,se$allele)]
gse <- se

idx <- 1:(ncol(gse)/2) * 2
cts <- assay(gse)
tot <- cts[,idx-1] + cts[,idx]
colnames(tot) <- sub("-a1","",colnames(gse)[idx])
plot(sort(colSums(tot)[22:112]))
# require 2e6 counts in the T-cells
cells_to_drop <- colnames(tot)[22:112][ colSums(tot)[22:112] < 2e6 ] 
alleles_to_drop <- sort(paste0(cells_to_drop,
                               rep(c("-a2","-a1"),each=length(cells_to_drop))))

library(fishpond)

gse <- gse[,!colnames(gse) %in% alleles_to_drop]
table(gse$cluster)/2
gse <- labelKeep(gse, minCount=1, minN=45)
table(mcols(gse)$keep)
gse <- swish(gse, x="allele", pair="cell", cov="cluster", interaction=TRUE)

hist(mcols(gse)$pvalue, breaks=40)
table(mcols(gse)$qvalue < .05)

suppressPackageStartupMessages(library(dplyr))
library(tibble)

res <- mcols(gse) %>%
  as.data.frame() %>%
  as_tibble() 

tot_sub <- tot[,!colnames(tot) %in% cells_to_drop]
res$mac <- rowSums(tot_sub[,1:21] >= 5)
res$tcell <- rowSums(tot_sub[,22:46] >= 5)
res$b6 <- rowSums(assay(gse)[,gse$allele == "a2"] >= 5)
res$cast <- rowSums(assay(gse)[,gse$allele == "a1"] >= 5)

res <- res %>%
  arrange(pvalue, -log10mean) %>%
  filter(mac >= 10 & tcell >= 10) %>%
  filter(b6 >= 5 & cast >= 5) %>%
  select(gene_id, gene_name, log10mean, log2FC, stat, pvalue, qvalue, 
         mac, tcell, b6, cast)

res_sub <- res %>%
  slice(1:20)

gse_symbol <- gse[res_sub$gene_id,]
rownames(gse_symbol) <- mcols(gse_symbol)$gene_name
gse_symbol <- computeInfRV(gse_symbol)
assays(gse_symbol) <- assays(gse_symbol)[!grepl("infRep",assayNames(gse_symbol))]

idx <- 1:(ncol(gse)/2) * 2
col_dat <- data.frame(condition=gse$cluster[idx],
                      row.names=sub("-a1","",colnames(gse)[idx]))

row_dat <- data.frame(absLFC=abs(res_sub$log2FC),
                      signLFC=factor(sign(res_sub$log2FC)),
                      minusLogQ=-log10(res_sub$qvalue),
                      row.names=res_sub$gene_name)

plotAllelicHeatmap(gse_symbol, idx=res_sub$gene_name,
                   annotation_row=row_dat,
                   annotation_col=col_dat,
                   cluster_rows=FALSE,
                   show_colnames=FALSE)

plotInfReps(gse_symbol, "Ly6e", x="allele", cov="cluster", thin=1)

tot_mat_to_plot <- log10(tot_sub[res_sub$gene_id,]+.1)
rownames(tot_mat_to_plot) <- res_sub$gene_name
pheatmap::pheatmap(tot_mat_to_plot,
                   annotation_col=col_dat,
                   cluster_rows=FALSE,
                   cluster_cols=FALSE,
                   show_colnames=FALSE,
                   main="Log10 total counts",
                   color=colorRampPalette(c("white","darkorchid1","darkorchid4"))(99))

# TSS level
load("allelic_counts_tss_filt.rda")
se$cluster <- factor(se$cluster)
levels(se$cluster) <- c("macrophage","T-cell")
se <- se[,order(se$cluster,se$cell,se$allele)]
ise <- se

ise <- ise[,!colnames(ise) %in% alleles_to_drop]
table(ise$cluster)/2
ise <- labelKeep(ise, minCount=1, minN=45)
table(mcols(ise)$keep)
ise <- swish(ise, x="allele", pair="cell", cov="cluster", interaction=TRUE)

tres <- mcols(ise) %>%
  as.data.frame() %>%
  as_tibble() 

tres_sub <- tres %>%
  dplyr::filter(gene_id %in% res_sub$gene_id) %>%
  dplyr::filter(qvalue < .5) %>%
  arrange(pvalue) %>%
  dplyr::select(gene_id, group_id, log2FC, pvalue, qvalue)

tres_sub %>%
  arrange(group_id) %>%
  inner_join(res[,c("gene_id","gene_name")])

suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
edb <- ah[["AH89211"]] # v102 EnsDb M.m.
ise_to_plot <- ise[,order(ise$allele, ise$cluster, ise$cell)]
ise_to_plot <- ise_to_plot[mcols(ise_to_plot)$keep,]
plotAllelicGene(ise_to_plot, gene="ENSMUSG00000022587", db=edb, cov="cluster")

plotAllelicGene(ise_to_plot, gene="ENSMUSG00000060803", db=edb, cov="cluster")
