# Bed file necessary columns 
# Chromosome ID [string] -> Chr
# Start genomic position of the phenotype -> Start 
# End genomic position of the phenotype -> End
# Phenotype ID -> ATAC_peak_1 Gene_id
# Phenotype group ID -> ATAC_peak_1 Gene_id
# Strand orientation [+/-] -> +

library(dplyr)
samples_list = read.table("metadata/sample_accession.txt")[,1]
#samples_list = c("ERS1199948", "ERS1199949", "ERS1199950", "ERS1199968")
filename = paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/featureCounts/", "ERS1199968", ".featureCounts.removed.txt", sep="")
counts = read.table(filename, header = TRUE)
counts_matrix = counts[,1:6]

for (samplename in samples_list){
  filename = paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/featureCounts/", samplename, ".featureCounts.removed.txt", sep="")
  counts = read.table(filename, header = TRUE)
  colnames(counts)[7] = "Fragmentcount"
  column = as.data.frame(counts$Fragmentcount)
  colnames(column) <- c(samplename)
  counts_matrix = dplyr::bind_cols(counts_matrix, column)
}

bed = dplyr::select(counts_matrix, dplyr::one_of("Chr", "Start", "End", "Geneid"))
bed = dplyr::mutate(bed, gid = Geneid)
bed = dplyr::bind_cols(bed, as.data.frame(counts_matrix[,5]))
colnames(bed)[1:6] <- c("#Chr", "start", "end", "pid", "gid", "strand")

bed = bind_cols(bed, counts_matrix[,7:length(colnames(counts_matrix))])

write.table(bed, "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/counts/counts_matrix.bed", sep="\t", row.names = FALSE, quote = FALSE )
