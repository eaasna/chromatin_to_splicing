#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install.packages("devtools")
library("devtools")
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")


filename = "counts/counts_matrix.bed"
#bed failis on columnames rea ees # see tuleb eemaldada
counts = read.table(filename)

accession = read.table("sample_name_genotype_full.txt")

colnames(counts)[1:6] <- c("Chr", "start", "end", "pid", "gid", "strand")
colnames(counts)[7:length(colnames(counts))] <- as.vector(accession[,2])

#gc = read.table(args[2])
gc = read.table("QC_measures/ATAC_GC_content.txt")
colnames(gc) <- c("Chr", "Start", "end", "Name", "Score", "Strand", "AT", "GC", "A", "C", "G", "T", "N", "Other", "Length")

library(dplyr)
metadata = dplyr::left_join(counts, gc, by="end") %>% dplyr::select(one_of("gid", "GC", "Length"))
colnames(metadata) <- c("gene_id", "percentage_gc_content", "length")
rownames(metadata) <- counts$gene_id

counts_matrix = data.matrix(counts[,7:length(colnames(counts))])
rownames(counts_matrix) <- counts$pid

#seqUtils meetod
library(cqn)
CQN_matrix = calculateCQN(counts_matrix, metadata)

meta = dplyr::left_join(counts, gc, by="end")[,1:6]

c = cbind(row.names(CQN_matrix), CQN_matrix)
colnames(c)[1] = "pid"

#add normalized counts, arranged by peak 
meta = full_join(meta, as.data.frame(c), by="pid")
colnames(meta)[1] <- "#Chr"

#outfile = args[3]
outfile = "normalized/counts_matrix_cqn.bed"
write.table(meta, outfile, sep="\t", row.names = FALSE, quote = FALSE)




