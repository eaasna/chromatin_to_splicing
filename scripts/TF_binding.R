library(readr)
library(Rsamtools)

# avatud kromatiin
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")

# r2 LD threshold
r = "08"

# transkriptsioonifaktori seondumissaidid
# ENCODE
# CTCF
#tf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz"
# HP1
#tf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/HP1/ENCFF417SVR.sorted.bed.gz" 
#tf = read.table(tf_file)
#bindingsites <- GRanges(seqnames = tf$V1, strand = c("*"), ranges = IRanges(start = tf$V2, end = tf$V3))
  
# QTLtools
# CTCF
tf_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz" # kõik seondumiskohad
tf = read.table(tf_file)
bindingsites <- GRanges(seqnames = tf$V2, strand = c("*"), ranges = IRanges(start = tf$V3, end = tf$V4, names = tf$V1))


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
  print(paste("H0: TF seondumine ei ole põhiline transkriptsiooni kontrollimehhanism"))
  print(paste("H1:", type, "seotud kühmude hulgas on TF seondumiskohad üleesindatud"))
  print(fisher.test(rbind(c(tabel[type, "eB"], tabel[type, "EB"]),c(tabel[type, "eb"],tabel[type, "Eb"])), alternative="less")$p.value)
}


#kas CTCF seondumissait kühmu juures tähendab, et kühm mõjutab suurema tõenäosusega fenotüüpi?
test_rida("featureCounts")
test_rida("upstream")
test_rida("contained")
test_rida("downstream")
test_rida("kokku")




# kas CTCFiga ülekattes olev caQTL mõjutab suurema tõenäosusega näiteks splaissimist kui geeniekspressiooni? 
test_tabel = function(tabel){
  #H0: mõlemat tüüpi RNA taseme tunnusega seotud kühmude hulgas leidub võrdne osakaal TF seondumissaite
  print("kas promootorites suurem osakaal kui featureCounts")
  print(fisher.test(rbind(c(tabel["featureCounts", "EB"], tabel["upstream", "EB"]),c(tabel["featureCounts", "Eb"], tabel["upstream", "Eb"])), alternative="less")$p.value)
  
  #H1: ekspressiooni ja 3' otsa kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
  print("kas featureCounts suurem osakaal kui promootoritel")
  print(fisher.test(rbind(c(tabel["upstream", "EB"], tabel["featureCounts", "EB"]),c(tabel["upstream", "Eb"], tabel["featureCounts", "Eb"])), alternative="less")$p.value)
  
  #H1: ekspressiooni ja splaissingu kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
  print("kas splaissimisel suurem osakaal kui featureCounts")
  print(fisher.test(rbind(c(tabel["featureCounts", "EB"], tabel["contained", "EB"]),c(tabel["featureCounts", "Eb"], tabel["contained", "Eb"])), alternative="less")$p.value)
  
  #H1: promootori ja splaissimise kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
  print("kas featureCounts suurem osakaal kui splaissimisel")
  print(fisher.test(rbind(c(tabel["contained", "EB"], tabel["featureCounts", "EB"]),c(tabel["contained", "Eb"], tabel["featureCounts", "Eb"])), alternative="less")$p.value)
  
  #H1: ekspressiooni ja splaissingu kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
  print("kas splaissimisel suurem osakaal kui promootoritel")
  print(fisher.test(rbind(c(tabel["upstream", "EB"], tabel["contained", "EB"]),c(tabel["upstream", "Eb"], tabel["contained", "Eb"])), alternative="less")$p.value)
  
  #H1: promootori ja splaissimise kühmude hulgas leidub erineval osakaalul CTCF seondumissaite
  print("kas promootoritel suurem osakaal kui splaissimisel")
  print(fisher.test(rbind(c(tabel["contained", "EB"], tabel["upstream", "EB"]),c(tabel["contained", "Eb"], tabel["upstream", "Eb"])), alternative="less")$p.value)
}

test_tabel(tabel)


