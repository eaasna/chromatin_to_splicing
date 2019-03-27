#install.packages("devtools")
library("devtools")
library(dplyr)
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

filename = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/counts/counts_matrix_PhaseIII.bed"
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
outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_fpkm.bed"

#remove -Inf values
mat = meta[,8:dim(meta)[2]]
mat[mat < -100] <- NA
meta[,8:dim(meta)[2]]=mat
complete = meta[complete.cases(meta),]

write.table(complete, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_fpkm_complete.bed", sep="\t", row.names = FALSE, quote = FALSE)
