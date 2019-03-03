# Rsquared.R loodud fail, mis näitab arvatavasti ühise lead variandiga piikide ja geenide paare

# kõik avatuse kühmud, millel on caQTL
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")

#r = '08'
r = '09'
#sub = 'featureCounts'
#sub = 'upstream'
#sub = 'contained'
sub = 'downstream'

# seotud avatuse kühmude ja geenide paarid
rsq_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r ,"/cqn_", sub, "_rsq.txt", sep="")
pairs = read.table(rsq_file)
colnames(pairs) <- c("peak_snp_id", "peak_snp_chr", "peak_snp_location", "gene_snp_id", "gene_snp_chr", "gene_snp_location", "rsq")

# kõik avatuse kühmud, millel on caQTL ja eQTL
enriched = dplyr::filter(caQTL, caQTL$V8 %in% pairs$peak_snp_id)
write.table(enriched, paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_", sub, "_enriched.txt", sep=""), col.names = F, quote = F, row.names = F, append = F, sep = '\t')


#Pärast seda ctcf.R