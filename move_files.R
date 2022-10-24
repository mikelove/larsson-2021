x <- read.delim("E-MTAB-12181.sdrf.txt")
names(x)[ncol(x)] <- "cell"
tab <- table(x$cell)
table(tab)
all(1:384 %in% names(tab))

for (i in 1:384) {
  print(i)
  x_sub <- x[x$cell == i,]
  cell_dir <- paste0("cell",i)
  system(paste("mkdir", paste0("fastq/",cell_dir)))
  files <- paste0(x_sub$Comment.ENA_RUN., ".fastq.gz")
  stopifnot(all(file.exists(files)))
  system(paste("mv",
               paste(files, collapse=" "),
               "-t",
               paste0("fastq/",cell_dir)))
}
