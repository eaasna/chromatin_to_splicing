#box plot of significant SNP with comparisons between genotypes
#jagan samplid selle järgi, kas genotype 0, 1, 2?
#normalized/counts_matrix_cqn.bed.gz <- siit võtan 
#genotypes/Kumasaka_100_samples.merged.vcf.gz

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
biocLite("GWASTools")
library("SNPRelate")
library("GWASTools")
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

#Convert VCF to GDS
snpgdsVCF2GDS("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.vcf.gz", 
              "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.gds", method="biallelic.only")
genotypes = gdsToMatrix("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.gds")
saveRDS(genotypes, "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.rds")


# Read 
genotypes = readRDS("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.rds")
atac = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_cqn.bed.gz", header = T)
QTLs = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL/FDR_corrected/cqn_permutations_100000_significant.txt", header = F)



# 1. vali välja mingi SNP + peak number
# 2. leia genotypes failist millistel samplitel selle SNP juures on 0, 1, 2
# 3. tee joonis nende gruppide kaupa atac jaotumisest

# QTLs with smallest p-value
most_significant = dplyr::arrange(QTLs, QTLs$V19)[1:5,]

#function to extract model p-value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

library(ggplot2)
plot_QTL <- function(SNP, peak, pos, ref, alt){
  library(data.table)
  smallest_genotype = setDT(as.data.frame(genotypes$genotypes[SNP,]), keep.rownames = TRUE)[]
  colnames(smallest_genotype) <- c("sample", "genotype")
  smallest_openness = setDT(as.data.frame(t(atac[which(atac$pid == peak), 7:length(colnames(atac))])), keep.rownames = TRUE)[]
  colnames(smallest_openness) <- c("sample", "atac")
  smallest = na.omit(dplyr::right_join(smallest_genotype, smallest_openness, by = "sample"))
  
  #write.table(smallest, paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL/plots/", SNP, ".txt"), row.names = F, quote = F)
  
  #lineaarne regressioon kromatiini avatuse ja genotüübi vahel
  fit = lm(formula = smallest$atac ~ smallest$genotype, data = smallest)
  
  plot = ggplot(data=smallest) + 
    geom_boxplot(mapping = aes(x = as.factor(genotype), y = atac, group = genotype)) + 
    #labs(x = SNP, y = "kromatiini avatus", title = paste("p-väärtus: ",lmp(fit), sep="")) + 
    labs(x = SNP, y = "kromatiini avatus", title = peak) + 
    scale_x_discrete(labels=c(paste(ref, ref, sep=""), paste(ref, alt, sep=""), paste(alt, alt, sep="")))
  #print(plot)
  ggsave(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL/plots/", SNP ,".png", sep=""), width = 5, height = 5)
}


SNP = "rs11242436"
peak = "ATAC_peak_223599"
# zgrep "rs11242436" Kumasaka_100_samples.merged.vcf.gz
# leia manuaalselt
ref = "A"
alt = "G"
plot_QTL(SNP, peak, pos, ref, alt)
# QTL genotüübist sõltuv kromatiini avatus ühe kühmu juures



