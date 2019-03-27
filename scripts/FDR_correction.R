#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# FDR correction teha ainult permutations failidel
sub = args[1]


sub = "downstream"
#tee ka contained ja downstream
#siis tagasi snakemake coloc
#full = read.table(paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/GEUVADIS/txrevise/naive.permuted.", sub ,".txt", sep = ""), header = F, stringsAsFactors = F)


full = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")

length(full$V19)

#grouped permutation run <- kasuta V21 so GEUVADIS
#vanilla permutation run <- kasuta V19
full=full[!is.na(full$V19),]
full$V19=p.adjust(full$V19, method = "fdr")

#write.table(full[which(full$V19 <= 0.1), ], paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/GEUVADIS/txrevise/naive.permuted.", sub,".significant.txt", sep = ""), quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(full[which(full$V19 <= 0.1), ], "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.txt.gz", quote=FALSE, row.names=FALSE, col.names=FALSE)


#cat /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL_rerun/cqn_permutations_*_100_100000_pca_5.txt >  /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL_rerun/cqn_100_100000_5_full.txt

perm = 100
windows = c("10000", "100000")
pcas = c(5, 10, 20)
for (window in windows){
  for (pca in pcas){
    print("number of PCA components")
    print(pca)
    print("window width")
    print(window)
    full = read.table(paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QTL_rerun/cqn_",perm,"_",window,"_",pca,"_full.txt", sep=""))
    full=full[!is.na(full$V19),]
    full$V19=p.adjust(full$V19, method = "fdr")
    print("number of significant QTLs")
    print(dim(full)[1])
    write.table(full[which(full$V19 <= 0.1), ], paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL_rerun/cqn_",perm,"_",window,"_",pca,"_significant.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


