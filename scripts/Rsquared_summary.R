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
with_peak08 = c(1164, 198, 344, 269)
with_peak09 = c(914, 136, 245, 178)
df = data.frame(total, with_peak08, with_peak09)
row.names(df) = c("featureCounts", "upstream", "contained", "downstream")



df = dplyr::mutate(df, without_peak = df$total - df$with_peak)
row.names(df) = c('featureCounts', 'upstream', 'contained', 'downstream')

write.table(df, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/Rsquared_summary.txt")




threshold = "with_peak09"
# H0: igat tüüpi eQTL'del on võrdselt ülekatet(coverage)
ttest = function(threshold, less, more) {
  print(paste("H1:", more, " rohkem ülekatet kui ", less))
  fisher.test(rbind(c(df[less, threshold], df[more, threshold]),c(df[less, "total"]-df[less, threshold],df[more, "total"]-df[more, threshold])), alternative="less")$p.value
}

ttest(threshold, "downstream", "upstream")
ttest(threshold, "contained", "upstream")
ttest(threshold, "featureCounts", "upstream")
ttest(threshold, "upstream", "featureCounts")
ttest(threshold, "contained", "featureCounts")
ttest(threshold, "downstream", "featureCounts")
ttest(threshold, "contained", "downstream")

