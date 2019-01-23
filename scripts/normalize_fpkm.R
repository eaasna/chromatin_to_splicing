#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install.packages("devtools")
library("devtools")
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

#filename = args[1]
filename = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/counts/counts_matrix.bed"
#bed failis on columnames rea ees # see tuleb eemaldada
counts = read.table(filename, header = TRUE)
counts_matrix = data.matrix(counts[,7:length(colnames(counts))])
rownames(counts_matrix) <- counts$pid

lengths = as.data.frame(counts$pid)
lengths = dplyr::bind_cols(lengths, dplyr::transmute(counts, length = end - start))
colnames(lengths) <- c("gene_id", "length")


#FPKM_matrix = calculateFPKM(counts_matrix, lengths)
# Siin on nimetus FPKM _matrix, tegelikult normaliseeritud TPM meetodil, mis on erinev!!! 
FPKM_matrix = calculateTPM(counts_matrix, lengths)
FPKM_matrix = log2(FPKM_matrix) + 1

meta = counts[,1:6]
meta = dplyr::bind_cols(meta, as.data.frame(FPKM_matrix))
colnames(meta)[1] <- "#Chr"
outfile = "normalized/counts_matrix_fpkm.bed"
accession = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QC_measures/sample_name_genotype_full.txt")

colnames(meta)[7:length(colnames(counts))] <- as.vector(accession[,2])
#outfile = args[2]
write.table(meta, outfile, sep="\t", row.names = FALSE, quote = FALSE)
