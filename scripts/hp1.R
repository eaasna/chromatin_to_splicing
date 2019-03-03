# r2 LD threshold
r = "09"

# transkriptsioonifaktori seondumissaidid
# ENCODE
# CTCF
#tf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz"
# HP1
tf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/HP1/ENCFF417SVR.sorted.bed.gz" 
tf = read.table(tf_file)
bindingsites <- GRanges(seqnames = tf$V1, strand = c("*"), ranges = IRanges(start = tf$V2, end = tf$V3))
  
# QTLtools
# CTCF
#tf_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz" # kõik seondumiskohad
#tf = read.table(tf_file)
#bindingsites <- GRanges(seqnames = tf$V2, strand = c("*"), ranges = IRanges(start = tf$V3, end = tf$V4, names = tf$V1))


RNA_tunnuse_rida = function(enriched, unenriched){
  unenriched_peaks = GRanges(seqnames = unenriched$V2, strand = c("*"), ranges = IRanges(start = unenriched$V3, end = unenriched$V4, names = unenriched$V1))
  unenriched_peaks_with_binding_site = subsetByOverlaps(unenriched_peaks, bindingsites)
  
  enriched_peaks = GRanges(seqnames = enriched$V2, strand = c("*"), ranges = IRanges(start = enriched$V3, end = enriched$V4, names = enriched$V1))
  enriched_peaks_with_binding_site = subsetByOverlaps(enriched_peaks, bindingsites)
  
  return(c(length(enriched_peaks_with_binding_site), length(enriched_peaks) - length(enriched_peaks_with_binding_site), length(unenriched_peaks_with_binding_site), length(unenriched_peaks) - length(unenriched_peaks_with_binding_site)))
}

tabel <- data.frame(matrix(ncol = 4, nrow = 5))
colnames(tabel) <- c("EB", "Eb", "eB", "eb")
row.names(tabel) <- c("featureCounts", "upstream", "contained", "downstream", "kokku")


# avatud kromatiin
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")


# RNA tunnused
featureCounts = read.table(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_featureCounts_enriched.txt", sep=""))
upstream = read.table(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_upstream_enriched.txt", sep=""))
contained = read.table(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r,"/cqn_contained_enriched.txt", sep=""))
downstream = read.table(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_downstream_enriched.txt", sep=""))

tabel["featureCounts", ] = RNA_tunnuse_rida(featureCounts, caQTL[which(!(caQTL$V1 %in% featureCounts$V1)),])
tabel["upstream", ] = RNA_tunnuse_rida(upstream, caQTL[which(!(caQTL$V1 %in% upstream$V1)),])
tabel["contained", ] = RNA_tunnuse_rida(contained, caQTL[which(!(caQTL$V1 %in% contained$V1)),])
tabel["downstream", ] = RNA_tunnuse_rida(downstream, caQTL[which(!(caQTL$V1 %in% downstream$V1)),])

kokku_unenriched = caQTL[which(!(caQTL$V1 %in% featureCounts$V1 | caQTL$V1 %in% upstream$V1 | caQTL$V1 %in% contained$V1 | caQTL$V1 %in% downstream$V1)),]
kokku_enriched = caQTL[which(caQTL$V1 %in% featureCounts$V1 | caQTL$V1 %in% upstream$V1 | caQTL$V1 %in% contained$V1 | caQTL$V1 %in% downstream$V1),]

tabel["kokku", ] = RNA_tunnuse_rida(kokku_enriched, kokku_unenriched)

print(tabel)



# kas seos RNA tunnusega mõjutab TF seondumiskohtade osakaalu?
test_rida = function(type){
  #print(paste("H0:", type, "seotud kühmude ja mitteseotud kühmude hulgas on võrdne osakaal TF seondumiskohti"))
  #print(paste("H1:", type, "seotud kühmude hulgas on rohkem TF seondumiskohti"))
  #print(fisher.test(rbind(c(tabel[type, "eB"],tabel[type, "EB"]),c(tabel[type, "eb"],tabel[type, "Eb"])), alternative="less")$p.value)

  print(paste("H0: TF seondumine ei ole põhiline transkriptsiooni kontrollimehhanism"))
  print(paste("H1: kühmud, mis sisaldavad TF seondumiskohta on suurema tõenäosusega seotud", type, "RNA QTLga"))
  print(fisher.test(rbind(c(tabel[type, "Eb"], tabel[type, "EB"]),c(tabel[type, "eb"],tabel[type, "eB"])), alternative="less")$p.value)
}

test_rida("featureCounts")
test_rida("upstream")
test_rida("contained")
test_rida("downstream")
test_rida("kokku")


# Enne võrdlesin ühte tüüpi RNA QTL'de puhul kas, leian CTCF mõju
# Nüüd võrdlen omavahel erinevat tüüpi RNA QTL'e
# andmed Tabel 2 aruandest
#H0: 
#H1: ekspressiooni ja promootori kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(227, 56),c(1170, 245)), alternative="less")$p.value

#H1: ekspressiooni ja splaissingu kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(227,76),c(1170,377)), alternative="less")$p.value

#H1: ekspressiooni ja 3' otsa kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(227,58),c(1170,283)), alternative="less")$p.value

#H1: promootori ja splaissimise kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(56,76),c(245,377)), alternative="less")$p.value

#H1: promootori ja 3' otsa kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(56,58),c(245,283)), alternative="less")$p.value

#H1: splaissimise ja 3' otsa kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
fisher.test(rbind(c(76,58),c(377,283)), alternative="less")$p.value
