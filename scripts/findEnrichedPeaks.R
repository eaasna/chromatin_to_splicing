# Rsquared.R loodud fail, mis näitab arvatavasti ühise lead variandiga piikide ja geenide paare
coverage = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_featureCounts_rsq.txt")
colnames(coverage) <- c("peak_snp_id", "peak_snp_chr", "peak_snp_location", "gene_snp_id", "gene_snp_chr", "gene_snp_location", "rsq")


# kõik piigid millel on caQTL
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")

# kõik piigid millel on caQTL ja eQTL
#enriched = dplyr::filter(caQTL, caQTL$V8 %in% coverage$peak_snp_id)

n = length(enriched[,1])

write.table(enriched, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_enriched.txt", col.names = F, quote = F, row.names = F, append = F, sep = '\t')

#Järgmisena Snakefile_ctcf, et tabix indekseerida z- pole vaja
#Pärast seda ctcf.R