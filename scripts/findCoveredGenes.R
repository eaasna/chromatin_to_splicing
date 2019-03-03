
pair_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_txrevise_downstream_rsq.txt"
#pair_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_txrevise_contained_rsq.txt"
#pair_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_txrevise_upstream_rsq.txt"
#pair_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_featureCounts_coverage.txt"

# Rsquared.R loodud fail, mis näitab arvatavasti ühise lead variandiga piikide ja geenide paare
pairs = read.table(pair_file)
colnames(pairs) <- c("peak_snp_id", "peak_snp_chr", "peak_snp_location", "gene_snp_id", "gene_snp_chr", "gene_snp_location", "rsq")


eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.downstream.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.contained.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.upstream.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/featureCounts.significant.sorted.txt.gz"

# kõik geenid, millel eQTL
eQTL = read.table(eqtl_file)

# kõik geenid mis on seotud ka kromatiini avatusega
coverage = dplyr::filter(eQTL, eQTL$V10 %in% pairs$gene_snp_id)

write.table(coverage, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_txrevise_downstream_coverage.txt", col.names = F, quote = F, row.names = F, append = F, sep = '\t')

#Järgmisena Snakefile_ctcf, et tabix indekseerida z- pole vaja
#Pärast seda ctcf.R