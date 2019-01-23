library(readr)
library(Rsamtools)


# ENCODE
#ctcf_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz"
#ctcf = read.table(ctcf_file)
#bindingsites <- GRanges(seqnames = ctcf$V1, strand = c("*"), ranges = IRanges(start = ctcf$V2, end = ctcf$V3))

# QTLtools
ctcf_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.txt.gz"
ctcf = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.txt.gz")
bindingsites <- GRanges(seqnames = ctcf$V2, strand = c("*"), ranges = IRanges(start = ctcf$V3, end = ctcf$V4, names = ctcf$V1))


# promootor
peak_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_upstream_enriched.sorted.txt"
# splaissing
#peak_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_contained_enriched.sorted.txt"
# 3 
# peak_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_downstream_enriched.sorted.txt"
# expression
#peak_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_featureCounts_enriched.sorted.txt"

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

# millised CTCF QTLd kattuvad teist sorti QTLdega
expression = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_upstream_enriched.txt")
ctcf = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_enriched.txt")
