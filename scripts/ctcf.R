library(readr)
library(Rsamtools)


# ENCODE
#ctcf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz"
#ctcf = read.table(ctcf_file)
#bindingsites <- GRanges(seqnames = ctcf$V1, strand = c("*"), ranges = IRanges(start = ctcf$V2, end = ctcf$V3))

# QTLtools
#ctcf_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.txt.gz" # seondumiskohad, millel leidub statistiliselt oluline QTL
ctcf_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz" # kõik seondumiskohad
ctcf = read.table(ctcf_file)
bindingsites <- GRanges(seqnames = ctcf$V2, strand = c("*"), ranges = IRanges(start = ctcf$V3, end = ctcf$V4, names = ctcf$V1))


#sub = "featureCounts"
#sub = "upstream"
#sub = "contained"
sub = "downstream"
r = "08"
# kühmud, mis on seotud teatud tüüpi RNA tunnusega
peak_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared", r, "/cqn_", sub ,"_enriched.txt", sep="")

# mul on vaja võtta ctcf seondumissaitide fail
# vaatan enriched kühme ja vaatan kui paljud neist on ctcf seondumissaidid
enriched = read.table(peak_file)
enriched_peaks <- GRanges(seqnames = enriched$V2, strand = c("*"), ranges = IRanges(start = enriched$V3, end = enriched$V4, names = enriched$V1))

enriched_peaks_with_binding_site = subsetByOverlaps(enriched_peaks, bindingsites)

# ATAC QTLtools väljund
caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")
caQTL = caQTL[which(!(caQTL$V1 %in% enriched$V1)),]

unenriched_peaks = GRanges(seqnames = caQTL$V2, strand = c("*"), ranges = IRanges(start = caQTL$V3, end = caQTL$V4, names = caQTL$V1))

# kui palju ülekatet kõigi ATAC kühmude ja CTCF kühmude vahel
unenriched_peaks_with_binding_site = subsetByOverlaps(unenriched_peaks, bindingsites)


# enriched + binding site +
EB = length(enriched_peaks_with_binding_site)
# enriched + binding site -
Eb = length(enriched_peaks) - length(enriched_peaks_with_binding_site)
# enriched - binding site +
eB = length(unenriched_peaks_with_binding_site)
# enriched - binding site -
eb = length(unenriched_peaks) - length(unenriched_peaks_with_binding_site)


# H0: enriched ja unenriched kühmude hulgas on võrdselt CTCF seondumiskohti
# H1: unenriched kühmude hulgas on rohkem CTCF seondumissaite
# unenriched on rohkem CTCF seondumiskohti?
fisher.test(rbind(c(EB, eB),c(Eb, eb)), alternative="less")$p.value

 
# H0: CTCF seondumiskoha esinemine kühmus ei mõjuta kas kühmul on eQTL
# H1: kühmud, mis sisaldavad CTCF seondumiskohta, mõjutavad geeniekspressiooni väiksema tõenäosusega

# CTCF seondumiskohaga kühmud mõjutavad ekspressiooni suurema tõenäosusega(on suurema tõenäosusega enriched)
fisher.test(rbind(c(Eb, EB),c(eb, eB)), alternative="less")$p.value

print(EB)
print(Eb)
print(eB)
print(eb)

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



# Nüüd võrdlen kõiki caQTL kühme mis on CTCF seondumissaidid nendega mis pole CTCF seondumissaidi
# Kas leidub erinevus enrichmentis üle kõigi RNA tunnuste?

caQTL = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")
all_peaks <- GRanges(seqnames = caQTL$V2, strand = c("*"), ranges = IRanges(start = caQTL$V3, end = caQTL$V4, names = caQTL$V1))
all_peaks_with_binding_site = subsetByOverlaps(all_peaks, bindingsites)


featureCounts = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_featureCounts_enriched.txt")
upstream = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_upstream_enriched.txt")
contained = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_contained_enriched.txt")
downstream = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared09/cqn_downstream_enriched.txt")

unenriched = caQTL[which(!(caQTL$V1 %in% featureCounts$V1 | caQTL$V1 %in% upstream$V1 | caQTL$V1 %in% contained$V1 | caQTL$V1 %in% downstream$V1)),]
unenriched_peaks = GRanges(seqnames = unenriched$V2, strand = c("*"), ranges = IRanges(start = unenriched$V3, end = unenriched$V4, names = unenriched$V1))

enriched = caQTL[which(caQTL$V1 %in% featureCounts$V1 | caQTL$V1 %in% upstream$V1 | caQTL$V1 %in% contained$V1 | caQTL$V1 %in% downstream$V1),]
enriched_peaks = GRanges(seqnames = enriched$V2, strand = c("*"), ranges = IRanges(start = enriched$V3, end = enriched$V4, names = enriched$V1))

enriched_peaks_with_binding_site = subsetByOverlaps(enriched_peaks, bindingsites)
unenriched_peaks_with_binding_site = subsetByOverlaps(unenriched_peaks, bindingsites)

# enriched + binding site +
EB = length(enriched_peaks_with_binding_site)
# enriched + binding site -
Eb = length(enriched_peaks) - length(enriched_peaks_with_binding_site)
# enriched - binding site +
eB = length(unenriched_peaks_with_binding_site)
# enriched - binding site -
eb = length(unenriched_peaks) - length(unenriched_peaks_with_binding_site)



print(EB)
print(Eb)
print(eB)
print(eb)

# kas CTCF seondumiskoha leidumine tähendab seda, et kühm on suurema tõenäosusega seotud mõne RNA taseme tunnusega
fisher.test(rbind(c(Eb, EB),c(eb, eB)), alternative="less")$p.value


# H1: unenriched kühmude hulgas on rohkem CTCF seondumissaite
fisher.test(rbind(c(EB, eB),c(Eb, eb)), alternative="less")$p.value


