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
cells_to_drop <- colnames(tot)[22:112][ colSums(tot)[22:112] < 2e6 ] # require 2e6 counts in the T-cells
alleles_to_drop <- sort(paste0(cells_to_drop, rep(c("-a2","-a1"),each=length(cells_to_drop))))

library(fishpond)

gse <- gse[,-which(colnames(gse) %in% alleles_to_drop)]
table(gse$cluster)/2
gse <- swish(gse, x="allele", pair="cell", cov="cluster", interaction=TRUE)

hist(mcols(gse)$pvalue, breaks=40)
table(mcols(gse)$qvalue < .05)

library(dplyr)
library(tibble)

res <- mcols(gse) %>%
  as.data.frame() %>%
  as_tibble() 

tot_sub <- tot[,-which(colnames(tot) %in% cells_to_drop)]
res$log10TotMac <- log10(rowSums(tot_sub[,1:21])+1)
res$log10TotTcell <- log10(rowSums(tot_sub[,22:46])+1)

res <- res %>%
  arrange(pvalue, -log10mean) %>%
  filter(log10TotMac > 3 & log10TotTcell > 3 & abs(log2FC) > 2) %>%
  select(gene_id, gene_name, log10TotMac, log10TotTcell,
         log10mean, log2FC, stat, pvalue, qvalue)

res_sub <- res %>%
  slice(1:20) %>%
  filter(!gene_name %in% c("Ptpra","Gpx1","Gm47428","Gm17494","Prex1","Cep85"))

gse_symbol <- gse[res_sub$gene_id,]
rownames(gse_symbol) <- mcols(gse_symbol)$gene_name

idx <- 1:(ncol(gse)/2) * 2
col_dat <- data.frame(condition=gse$cluster[idx],
                      row.names=sub("-a1","",colnames(gse)[idx]))

row_dat <- data.frame(absLFC=abs(res_sub$log2FC),
                      signLFC=factor(sign(res_sub$log2FC)),
                      log10TotTcell=res_sub$log10TotTcell,
                      log10TotMac=res_sub$log10TotMac,
                      row.names=res_sub$gene_name)

plotAllelicHeatmap(gse_symbol, idx=res_sub$gene_name,
                   annotation_row=row_dat,
                   annotation_col=col_dat,
                   cluster_rows=FALSE,
                   show_colnames=FALSE)

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

ise <- ise[,-which(colnames(ise) %in% alleles_to_drop)]
table(ise$cluster)/2
ise <- swish(ise, x="allele", pair="cell", cov="cluster", interaction=TRUE)

tres <- mcols(ise) %>%
  as.data.frame() %>%
  as_tibble() 

tres_sub <- tres %>%
  dplyr::filter(gene_id %in% res_sub$gene_id) %>%
  dplyr::filter(qvalue < .05) %>%
  arrange(pvalue) %>%
  dplyr::select(gene_id, group_id, log2FC, pvalue, qvalue)

library(AnnotationHub)
ah <- AnnotationHub()
edb <- ah[["AH89211"]] # v102 EnsDb M.m.
ise_to_plot <- ise[,order(ise$allele, ise$cluster, ise$cell)]
plotAllelicGene(ise_to_plot, gene="ENSMUSG00000022587", db=edb, cov="cluster")

