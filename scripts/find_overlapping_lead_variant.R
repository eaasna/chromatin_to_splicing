#permutation run ühiste lead variantide leidmine
expression = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_10000_significant.sorted.txt.gz", header = F)
featureCounts = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/featureCounts.significant.sorted.txt.gz", header = F)
#colnames(fpkm) <- c("phenotype_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "n_snps", "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "p_nominal", "beta", "is_lead")
#colnames(featureCounts) <- c("phenotype_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "n_snps", "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "p_nominal", "beta", "is_lead")

colnames(expression)[9:10] <- c('snp_chr', 'snp_start')
colnames(featureCounts)[11:12] <- c('snp_chr', 'snp_start')

overlap = dplyr::inner_join(expression, featureCounts, by = c("snp_chr" = "snp_chr", "snp_start" = "snp_start"))[,c(1, 2, 3, 4, 6, 8, 10, 19, 20, 21, 22, 23, 27, 29, 30, 38)]
colnames(overlap) <- c("peak_id", "peak_chr", "peak_start", "peak_end", "peak_n_snps", "peak_snp_id", "peak_snp_location", "peak_p_beta", "gene_id", "gene_chr", "gene_start", "gene_end", "gene_n_snps", "gene_snp_id", "gene_snp_location", "gene_p_beta")


write.table(overlap, file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/shared_lead/cqn.10000.featureCounts.shared.lead.txt", col.names = F, quote = F, row.names = F, append = F)

#kontroll
fpkm = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/shared_lead/fpkm.10000.featureCounts.shared.lead.txt", header=F)
cqn = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/shared_lead/cqn.10000.featureCounts.shared.lead.txt", header=F)
