source("https://bioconductor.org/biocLite.R")
library("SNPRelate")
library("GWASTools")
library("rtracklayer")
library("GenomicRanges")
library("devtools")
load_all("seqUtils")
#biocLite("SNPRelate")
#library("SNPRelate")

chrom = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")

genotypeList = list()

for (c in chrom){
  #Make a GDS file
  snpgdsVCF2GDS(paste("genotypes/Kumasaka_100_samples.chr", c, ".vcf.gz", sep=""),
                paste("genotypes/Kumasaka_100_samples.chr", c, ".gds", sep=""),
              method="biallelic.only") #method = "copy.num.of.ref"

  #Make a RDS file
  genotypes = gdsToMatrix(paste("genotypes/Kumasaka_100_samples.chr", c, ".gds", sep=""))
  saveRDS(genotypes, paste("genotypes/Kumasaka_100_samples.chr", c, ".rds", sep=""))
  
  #Siit saan SNP asukohad
  genotypeRDS = readRDS(paste("genotypes/Kumasaka_100_samples.chr", c, ".rds", sep=""))
  
  #Make GRanges object of SNP locations for each chromosome
  gr_SNP <- makeGRangesFromDataFrame(genotypeRDS$snpspos,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field="chr",
                                     start.field="pos",
                                     end.field="pos")
  
  genotypeList <- c(genotypeList, gr_SNP)
}

#GRanges object of all SNP positions
SNP_positions = do.call("c", genotypeList)

#Siit saan peak asukohad
gtf <- rtracklayer::import('annotations/ATAC_joined_peaks.gtf')
gtf_df=as.data.frame(gtf)

gtf_df = dplyr::mutate(gtf_df, start = gtf_df$start + round(gtf_df$width/2) - 50000, end = gtf_df$end - round(gtf_df$width/2) + 50000)
gtf_df$start[gtf_df$start < 0] <- 0
gtf_df <- dplyr::select(gtf_df, dplyr::one_of(c("seqnames", "start", "end", "width", "gene_id")))
colnames(gtf_df)[5] <- "peak_id"

#Make peak GRanges object
gr_peak <- makeGRangesFromDataFrame(gtf_df,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field="seqnames",
                                    start.field="start",
                                    end.field="end")

overlaps <- findOverlaps(gr_peak, SNP_positions)

peaks_with_overlap <- subsetByOverlaps(gr_peak, SNP_positions)$peak_id
all_peaks <- gtf_df$peak_id
peaks_witout_overlap <- setdiff(all_peaks, peaks_with_overlap)
write(peaks_witout_overlap, file="peaks_with_no_nearby_SNPs.txt")
