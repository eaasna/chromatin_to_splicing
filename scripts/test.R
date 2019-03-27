full = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QTL/cqn_permutations_100000_full.txt.gz", header = F, stringsAsFactors = F)

#png(filename=paste(path, "plots/", args[1] ,".png", sep=""))
sub = "ilus_cqn"

png(filename=paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL/plots/", sub, ".png", sep=""))
#plot(full$V18, full$V19, xlab="Direct method", ylab="Beta approximation", main=args[1])
#plot(full$V18, full$V19, xlab="Otsene meetod", ylab="Beta ennustus", main=sub)
#abline(0, 1, col="red")
hist(full$V19, xlab="p-väärtus", ylab="sagedus", breaks = 100, main="Permuteeritud p-väärtuste histogramm")
dev.off()

#file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_downstream_rsq.txt"
#file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_contained_rsq.txt"
#file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_upstream_rsq.txt"
file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_featureCounts_rsq.txt"

pairs = read.table(file)
pairs = dplyr::transmute(pairs, distance = abs(pairs$V3 - pairs$V6))

hist(pairs$distance, breaks = 200, main = "featureCounts")

# 
cqn = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/normalized/counts_matrix_cqn.bed.gz", header = TRUE)
fpkm = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/normalized/counts_matrix_fpkm.bed.gz", header = TRUE)

full = dplyr::left_join(fpkm, cqn, by = "pid")
cqn = full[1:106]
fpkm = full[,107:211]


apply(full, 1, function(x) print(cor(as.factor(x[,7:106]), as.factor(x[,112:211]))))


max(fpkm$HG00096.y)
min(fpkm$HG00096.y)

library(ggplot2)
plot = ggplot(data=full) + 
  geom_jitter(aes(x = cqn$HG00096.x, y = fpkm$HG00096.y)) + 
  #labs(x = SNP, y = "kromatiini avatus", title = paste("p-väärtus: ",lmp(fit), sep="")) + 
  labs(x = "CQN", y = "FPKM", title = "1 proov")
print(plot)


# Bed file necessary columns 
# Chromosome ID [string] -> Chr
# Start genomic position of the phenotype -> Start 
# End genomic position of the phenotype -> End
# Phenotype ID -> ATAC_peak_1 Gene_id
# Phenotype group ID -> ATAC_peak_1 Gene_id
# Strand orientation [+/-] -> +


ctcf = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")
ctcf_bed = data.frame(matrix(ncol = 4, nrow = dim(ctcf)[1]))
ctcf_bed = ctcf[,c(2,3, 4)]

write.table(ctcf_bed,"/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/CTCF/CTCF.bed", col.names = F, quote = F, row.names = F, append = F, sep = '\t')




