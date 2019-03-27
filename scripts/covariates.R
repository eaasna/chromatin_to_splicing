genotype = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.pca", header = TRUE)

norm_types = c("cqn")
for (norm in norm_types){
  normalization = read.table(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/",norm,".pca", sep=""), header = TRUE)  
  cov_counts = c(5, 10, 20)
  for (count in cov_counts){
    cov = dplyr::bind_rows(head(genotype, count), head(normalization, count))
    write.table(cov, paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/covariates/",norm,"_",count,".pca",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  }
}


#vaja on luua sellised covariates failid millel päises HG tüüpi identifikaatorid
