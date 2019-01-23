match = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/snpsnap/SNPsnap_upstream/matched_snp.txt", sep=' ')
#match = match[match$pos %in% genotype$snpspos$pos, ]
colnames(match) <- c('chr', 'pos')

genotype = readRDS("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.rds")
snppos = genotype$snpspos

match = dplyr::mutate(match, chr = as.character(chr))
match = dplyr::inner_join(match, snppos, by = c("pos", "chr"))


annotated = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/snpsnap/SNPsnap_upstream/matched_snps_annotated.txt", sep="\t", header = T)


#kõik ctcf seondumissaidid
feature = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.bed.gz")


