library("rtracklayer")
library("devtools")
library("dplyr")
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

#Load all sample names from disk
#sample_names = read.table("/gpfs/hpchome/evelin95/metadata/sample_name_genotype_phaseIII.txt")[,1]
sample_names = read.table("sample_name_genotype_phaseIII.txt")[,1]
#sample_names = c("ERS1199948", "ERS1199949", "ERS1199950")

#Import peak call
#Find peaks
peak_list = loadNarrowPeaks("peaks", sample_names)
peak_count = lapply(peak_list, length) %>%
  plyr::ldply(.id = "sample_id") %>%
  dplyr::rename(peak_count = V1)

#Make peak list
# saveRDS(peak_list, "results/ATAC/ATAC_peak_list.rds")
# peak_list = readRDS("results/ATAC/ATAC_peak_list.rds")

#Filter peaks by overlap
filtered = filterOverlaps(peak_list, 3)
all_peaks = purrr::reduce(filtered, union.Vector)

#Export as a gff file
# rtracklayer::export.gff3(all_peaks, "annotations/ATAC_consensus_peaks.gff3")
# rtracklayer::export.bed(all_peaks, "annotations/ATAC_consensus_peaks.bed")

# Join all peaks together
joined_peaks = GenomicRanges::reduce(all_peaks)

#Add additional columns
metadata = data_frame(type = "exon", gene_id = paste("ATAC_peak_", c(1:length(joined_peaks)), sep = ""))
elementMetadata(joined_peaks) = metadata
strand(joined_peaks) = "+"

rtracklayer::export.gff3(joined_peaks, "annotations/ATAC_joined_peaks.gff3")
rtracklayer::export.bed(joined_peaks, "annotations/ATAC_joined_peaks.bed")


#Use peak count per sample as QC measures
write.table(peak_count, "QC_measures/macs2_peaks_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(peak_list, "results/ATAC/ATAC_peak_list.rds")

#enne seda kopeeri uus loodud gff3 fail rtracklayer kausta
?import.gff3
test_path <- system.file("tests", package = "rtracklayer")
test_gff3 <- file.path(test_path, "ATAC_joined_peaks.gff3")
test <- import(test_gff3)
export(test,"ATAC_joined_peaks.gtf","gtf")
