clust_id <- read.table("Rmexbursting_ClusterIdentity.txt")
names(clust_id) <- c("fastq","cluster")
metadata <- read.delim("E-MTAB-12181.sdrf.txt")
metadata <- metadata[grep("P8467", metadata$Scan.Name),
                     c("Scan.Name","Factor.Value.single.cell.identifier.")]
names(metadata) <- c("fastq","cell")
metadata$fastq <- sub("(P8467_.*?)_.*","\\1",metadata$fastq)
metadata <- metadata[!duplicated(metadata$cell),]
metadata <- metadata[order(metadata$cell),]
metadata <- merge(metadata, clust_id)
write.csv(metadata, file="metadata.csv", row.names=FALSE, quote=FALSE)
