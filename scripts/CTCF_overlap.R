library(readr)
library(Rsamtools)
library(GenomicRanges)

encode_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz"
encode = read.table(encode_file)
encode_bindingsites <- GRanges(seqnames = encode$V1, strand = c("*"), ranges = IRanges(start = encode$V2, end = encode$V3))

yuri_file = "/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz" # kõik seondumiskohad
yuri = read.table(yuri_file)
yuri_bindingsites <- GRanges(seqnames = yuri$V2, strand = c("*"), ranges = IRanges(start = yuri$V3, end = yuri$V4, names = yuri$V1))

hits = findOverlaps(encode_bindingsites, yuri_bindingsites)
overlaps <- pintersect(encode_bindingsites[queryHits(hits)], yuri_bindingsites[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(yuri_bindingsites[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]


hits = findOverlaps(yuri_bindingsites, encode_bindingsites)
overlaps <- pintersect(yuri_bindingsites[queryHits(hits)], encode_bindingsites[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(encode_bindingsites[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]