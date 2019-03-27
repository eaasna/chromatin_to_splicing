
#install.packages("devtools")
library("devtools")
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

filename = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/counts/counts_matrix_PhaseIII.bed"
counts = read.table(filename, header = TRUE)

gc = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QC_measures/ATAC_GC_content.txt")
colnames(gc) <- c("Chr", "start", "end", "Name", "Score", "Strand", "AT", "GC", "A", "C", "G", "T", "N", "Other", "Length")
gc$start = gc$start + 1

library(dplyr)
metadata = dplyr::left_join(counts, gc, by=c("start"="start","end"="end")) %>% dplyr::select(one_of("gid", "GC", "Length"))
colnames(metadata) <- c("gene_id", "percentage_gc_content", "length")
rownames(metadata) <- counts$gene_id

counts_matrix = data.matrix(counts[,7:length(colnames(counts))])
rownames(counts_matrix) <- counts$pid

#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("cqn"))

library(cqn)
#seqUtils meetod
CQN_matrix = calculateCQN(counts_matrix, metadata)

#meta = dplyr::left_join(counts, gc, by=c("start"="start","end"="end"))[,1:6]
meta = counts[,1:6]

c = cbind(row.names(CQN_matrix), CQN_matrix)
colnames(c)[1] = "pid"

#add normalized counts, arranged by peak 
meta = full_join(meta, as.data.frame(c), by="pid")
colnames(meta)[1] <- "#Chr"

outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_cqn.bed"
write.table(meta, outfile, sep="\t", row.names = FALSE, quote = FALSE)




