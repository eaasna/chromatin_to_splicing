library(readr)
library(Rsamtools)

# siin kõikide regulatoorsete alade annotatsioonid 
# allikas http://genome.ucsc.edu/cgi-bin/hgFileUi?g=wgEncodeBroadHmm&db=hg18
# annotatsioonid sellest artiklist file:///C:/Users/evelin95/Zotero/storage/ZH3X37NV/Ernst%20et%20al.%20-%202011%20-%20Mapping%20and%20analysis%20of%20chromatin%20state%20dynamics%20i.pdf
annotations = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/annotations/wgEncodeBroadHmmGm12878HMM.sorted.bed.gz", sep = "\t")

# kõik ATAC kühmud millel leidub significant caQTL
peaks = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz")

# huvipakkuvad kühmud
#peaks = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_contained_enriched.sorted.txt.gz")
gr_peaks<- GRanges(seqnames = peaks$V2, strand = c("*"), ranges = IRanges(start = peaks$V3, end = peaks$V4, names = peaks$V1))


regulatory_df = data.frame(row.names = unique(annotations$V4), "peaks" = vector(mode = "numeric", length = 15), "CTCF_sites" = vector(mode = "numeric", length = 15))

#ENCODE
#ctcf = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/ENCFF960ZGP.sorted.bed.gz")
#gr_bindingsites <- GRanges(seqnames = ctcf$V1, strand = c("*"), ranges = IRanges(start = ctcf$V2, end = ctcf$V3))
#QTLtools 
ctcf = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")
gr_bindingsites <- GRanges(seqnames = ctcf$V2, strand = c("*"), ranges = IRanges(start = ctcf$V3, end = ctcf$V4, names = ctcf$V1))



# CTCF binding sites overlaping peaks
binding_sites_within_peaks = subsetByOverlaps(gr_bindingsites, gr_peaks)

peaks_within_binding_sites = subsetByOverlaps(gr_peaks, gr_bindingsites)
#kui suur osa kühmudest sisaldab ctcf seondumiskohta
length(gr_peaks)/length(unique(names(peaks_within_binding_sites)))

#huvipakkuvate kühmude arv
print(length(peaks$V1))

for (type in unique(annotations$V4)){
  annotation_subset = dplyr::filter(annotations, annotations$V4 == type)
  
  # nüüd leia kui paljud kühmud kuuluvad regulatoorsetesse aladesse
  gr_annotation <- GRanges(seqnames = annotation_subset$V1, strand = c("*"), ranges = IRanges(start = annotation_subset$V2, end=annotation_subset$V3))
  gr_peaks <- GRanges(seqnames = peaks$V2, strand = c("*"), ranges = IRanges(start = peaks$V3, end=peaks$V4, names = peaks$V1))
  
  #peaks that overlap with "type" region
  peak_subset <- subsetByOverlaps(gr_peaks, gr_annotation)
  
  regulatory_df[type,"peaks"] = length(names(ranges(peak_subset)))
  #write.table(names(ranges(peak_subset)), paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ATAC/cqn_contained_enriched_", type,".txt", sep = "")) 

  #CTCF binding sites overlaping peaks that overlap "type" regulatory region
  annotation_subset = subsetByOverlaps(gr_annotation, gr_peaks) # promoters/enhansers within peaks that have caQTL
  regulatory_df[type,"CTCF_sites"] = length(ranges(subsetByOverlaps(binding_sites_within_peaks, annotation_subset)))
}

regulatory_df

# huvipakkuvate kühmudega kattuvate CTCF seondumiskohtade arv
print(length(ranges(binding_sites_within_peaks)))


# mitu caQTL asuvad mõjutatava kühmu sees?
dim(peaks[which(peaks$V10 >= peaks$V3 & peaks$V10 <= peaks$V4), ])


hist(abs((wide$V4 - wide$V3)/2 + wide$V3 - wide$V10), breaks = 100, main = "", ylab="Sagedus", xlab = "QTL kaugus kühmu keskpunktist")


