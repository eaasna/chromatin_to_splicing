outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_downstream_coverage.txt"
#outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_contained_coverage.txt"
#outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_upstream_coverage.txt"
#outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_featureCounts_coverage.txt"

eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.downstream.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.contained.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.upstream.sorted.txt.gz"
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/featureCounts.significant.sorted.txt.gz"




# võta kogu geenide arv millel on eQTL
# võta geenide arv, millel on üks kühmuga ühine eQTL > 0.8

rsq = read.table(outfile, header =F )
qtl = read.table(eqtl_file, header = F)

nrow(rsq)
length(unique(qtl$V1))



total = c(5957, 1200, 2716, 2065)
with_peak = c(1164, 198, 344, 269)
df = data.frame(total, with_peak)


#sum = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/Rsquared_summary.txt", sep = ',', header = T)

df = dplyr::mutate(df, without_peak = df$total - df$with_peak)
row.names(df) = c('featureCounts', 'upstream', 'contained', 'downstream')

write.table(df, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/Rsquared_summary.txt")


# H0: igat tüüpi eQTL'del on võrdselt ülekatet(coverage)
# H1: promootoril(upstream) on rohkem ülekatet kui 3' otsal
fisher.test(rbind(c(269,198),c(1796,1002)), alternative="less")$p.value

# H1: promootoril on rohkem ülekatet kui keskmisel eksonil
fisher.test(rbind(c(344,198),c(2372,1002)), alternative="less")$p.value

# H1: promootoril on rohkem ülekatet kui featureCounts
fisher.test(rbind(c(1164,198),c(4793,1002)), alternative="less")$p.value

# H1: featureCounts on rohkem ülekatet kui promootoril
fisher.test(rbind(c(198,1164),c(1002,4793)), alternative="less")$p.value

# H1: featureCounts on rohkem ülekatet kui 3' otsal
fisher.test(rbind(c(269,1164),c(1796,4793)), alternative="less")$p.value

# H1: featureCounts on rohkem ülekatet kui keskmisel eksonil
fisher.test(rbind(c(344,1164),c(2372,4793)), alternative="less")$p.value

# H1: keskmisel eksonil on rohkem ülekatet kui 3' otsal
fisher.test(rbind(c(269,344),c(1796,2372)), alternative="less")$p.value


txrevise = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/GEUVADIS/txrevise/naive.permuted.upstream.significant.txt")
snap = unique(txrevise$V10)
write.table(snap, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/GEUVADIS/txrevise/upstream.significant.snps.txt", col.names = F, quote = F, row.names = F, append = F)

