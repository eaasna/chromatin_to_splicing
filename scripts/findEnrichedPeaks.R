# Rsquared.R loodud fail, mis näitab arvatavasti ühise lead variandiga piikide ja geenide paare

# kõik avatuse kühmud, millel on caQTL
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")


r = '08'
#r = '09'


find_enriched = function(rsq_file, outfile){
  pairs = read.table(rsq_file)
  colnames(pairs) <- c("peak_snp_id", "peak_snp_chr", "peak_snp_location", "gene_snp_id", "gene_snp_chr", "gene_snp_location", "rsq")
  # kõik avatuse kühmud, millel on caQTL ja eQTL
  enriched = dplyr::filter(caQTL, caQTL$V8 %in% pairs$peak_snp_id)
  write.table(enriched, outfile, col.names = F, quote = F, row.names = F, append = F, sep = '\t')
  
}

# kromatiini kühmud, mis on seotud RNA tunnusega
for (sub in c('featureCounts', 'upstream', 'contained', 'downstream')){
  rsq_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r ,"/cqn_", sub, "_rsq.txt", sep="")
  outfile = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_", sub, "_enriched.txt", sep="")
  find_enriched(rsq_file, outfile)
}


# kromatiini kühmud, mis on seotud mõne CTCF QTLga
rsq_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_rsq", r ,".txt", sep="")
outfile = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_enriched", r, ".txt", sep="")
find_enriched(rsq_file, outfile)


# kromatiini kühmud, mis on seotud RNA tunnusega ja CTCF QTLga
for (sub in c('featureCounts', 'upstream', 'contained', 'downstream')){
  rsq_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_", sub, "_rsq", r ,".txt", sep="")
  outfile = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_", sub, "_enriched", r, ".txt", sep="")
  find_enriched(rsq_file, outfile)
}





