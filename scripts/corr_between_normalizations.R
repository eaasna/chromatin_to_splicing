filename = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/counts/counts_matrix.bed"
#bed failis on columnames rea ees # see tuleb eemaldada
counts = read.table(filename, header = TRUE)
dim(counts)
not_phaseIII = c("ERS1200001", "ERS1200008", "ERS1200012", "ERS1200014", "ERS1200018", "ERS798631", "ERS798637", "ERS798626", "ERS798629")
counts = dplyr::select(counts, -not_phaseIII)

#filtering out 9 samples not in Phase III
accession = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QC_measures/run_sample_accession_PhaseIII.txt")
colnames(counts)[1:6] <- c("Chr", "start", "end", "pid", "gid", "strand")
colnames(counts)[7:length(colnames(counts))] <- as.vector(accession[,2])

write.table(counts, "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/counts/counts_matrix_PhaseIII.bed", sep="\t", row.names = FALSE, quote = FALSE)

cqn = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/normalized/counts_matrix_cqn.bed")
fpkm = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/normalized/counts_matrix_fpkm_complete.bed")

c = c("V25", "V35", "V45", "V55")

for (col in c){
  corr_df = dplyr::left_join(fpkm[,c("V4", col)], cqn[,c("V4", col)], by="V4")
  print(cor(corr_df[,2], corr_df[,3]))  
}
