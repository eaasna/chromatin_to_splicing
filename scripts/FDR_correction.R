#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# FDR correction teha ainult permutations failidel
sub = args[1]


sub = "downstream"
#tee ka contained ja downstream
#siis tagasi snakemake coloc
#full = read.table(paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/GEUVADIS/txrevise/naive.permuted.", sub ,".txt", sep = ""), header = F, stringsAsFactors = F)

full = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")


#grouped permutation run <- kasuta V21 so GEUVADIS
#vanilla permutation run <- kasuta V19
full=full[!is.na(full$V19),]
full$V19=p.adjust(full$V19, method = "fdr")

#write.table(full[which(full$V19 <= 0.1), ], paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/GEUVADIS/txrevise/naive.permuted.", sub,".significant.txt", sep = ""), quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(full[which(full$V19 <= 0.1), ], "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.txt.gz", quote=FALSE, row.names=FALSE, col.names=FALSE)
