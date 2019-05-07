# tegelen Yuri CTCF andmetega
# ainult 0.8 r2 väärtusi
# ainult contained ehk splaissimise QTLd


library("dplyr")
library("readr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library('Rsamtools')
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")


path = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/"

  

#Selles piirkonnas on kühmud ATAC_peak_257613(rs2429285), ATAC_peak_257614(rs1799127) seotud CTCF SNP'ga (rs573614910 pos 76515459). caQTL ise on väljaspool kujutatud piirkonda. Aga splaiss QTL????
region_coords = c(76514000,76518000)
rs_id = "rs73142003"
gene_id = 'chr7_76514308_G_A'
#vcf_file$snpspos[which(vcf_file$snpspos$chr==7 &vcf_file$snpspos$pos>76516418),]
#vcf_file$snpspos[which(vcf_file$snpspos$chr==7 & vcf_file$snpspos$pos>76514000 & 76518000 > vcf_file$snpspos$pos),]



region_coords = c(76461000,76506000)
vcf_file$snpspos[which(vcf_file$snpspos$chr==7 & vcf_file$snpspos$pos>region_coords[1] & region_coords[2] > vcf_file$snpspos$pos),]
gene_id = "chr7_76462814_C_T"



region_coords = c(76483500,76484000)
vcf_file$snpspos[which(vcf_file$snpspos$chr==7 & vcf_file$snpspos$pos>region_coords[1] & region_coords[2] > vcf_file$snpspos$pos),"snpid"][[1]]
gene_id = "chr7_76462814_C_T"

#----------------------------------------------------------------------------------------------#
#'
#'CTCF andmete sisse lugemine
#'
#

#kõik seondumiskohad
ctcf = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")
ctcf_ranges <- GRanges(seqnames = ctcf$V2, strand = c("*"), ranges = IRanges(start = ctcf$V3, end = ctcf$V4, names = ctcf$V1))

#r2 seotud splaiss ja CTCF lookuste paarid
ctcf_contained = read.table(paste(path,"results/ctcf/cqn_contained_rsq08.txt", sep=""))

#seondumiskohad, mille QTL on seotud splaissimisega
ctcf_contained_binding = ctcf[which(ctcf$V8 %in% ctcf_contained$V4),]
ctcf_contained_ranges = GRanges(seqnames = ctcf_contained_binding$V2, strand = c("*"), ranges = IRanges(start = ctcf_contained_binding$V3, end = ctcf_contained_binding$V4, names = ctcf_contained_binding$V1))

CTCF_peak_metadata = ctcf[,c("V1","V2", "V3","V4", "V5")]
colnames(CTCF_peak_metadata) = c("pid","chr","start", "end","strand")


CTCF_peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 7, CTCF_peak_metadata)
CTCF_peak_annot$peak_annot$transcript_id='CTCF'
CTCF_peak_annot$peak_annot$gene_id = 'CTCF'
CTCF_peak_annot$peak_annot$gene_name = 'CTCF CHiP-seq'
names(CTCF_peak_annot$peak_list) = 'CTCF'

#//////////////////////////
#geenide asukoht

if (FALSE){
  CTCF_peak_plot = plotTranscripts(CTCF_peak_annot$peak_list, CTCF_peak_annot$peak_list, CTCF_peak_annot$peak_annot, rescale_introns = FALSE, 
                                   region_coords = region_coords, connect_exons = FALSE, transcript_label = FALSE)
}

#//////////////////////////

rm(ctcf, ctcf_ranges)

#----------------------------------------------------------------------------------------------#
#'
#'RNA seq andmete sisse lugemine
#'

#proovid
RNA_sample_metadata = read.table("/gpfs/hpchome/a72094/hpc/datasets/controlled_access/SampleArcheology/studies/cleaned/GEUVADIS.tsv", header = TRUE)[,c("sample_id","genotype_id","condition")]
colnames(RNA_sample_metadata)[3] = "condition_name"


# RNA seq lugemite arv used for calculating library size
RNA_counts = read.table("/gpfs/hpchome/a72094/hpc/projects/RNAseq_pipeline/results/expression_matrices/featureCounts/GEUVADIS.tsv.gz", header = TRUE)

RNA_meta_df = wiggleplotrConstructMetadata(RNA_counts, RNA_sample_metadata, "/gpfs/hpc/home/a72094/projects/RNAseq_pipeline/processed/GEUVADIS/bigwig", bigWig_suffix = ".str1.bw", condition_name_levels = c("naive"))


# vcf to gds file
library(SNPRelate)
gds_file = paste0(path, "genotypes/chr7.gds")
snpgdsVCF2GDS(paste0(path, "genotypes/7.vcf"), gds_file)
vcf_file=gdsToMatrix(gds_file)


# RNA metadata
RNA_peak_metadata = read.table("/gpfs/hpchome/a72094/hpc/projects/RNAseq_pipeline/metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz", header = TRUE)[,c("phenotype_id", "chromosome","gene_start","gene_end")]
colnames(RNA_peak_metadata) = c("pid","chr","start", "end")
RNA_peak_metadata$strand = rep("*",dim(RNA_peak_metadata)[1])

#/////////////////
if (FALSE){
  #Fetch all peaks in the region and make peak annot plot
  RNA_peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 7, RNA_peak_metadata)
  RNA_peak_annot$peak_annot$transcript_id='RNA'
  RNA_peak_annot$peak_annot$gene_id = 'RNA'
  RNA_peak_annot$peak_annot$gene_name = 'RNA-seq'
  names(RNA_peak_annot$peak_list) = 'RNA'
  
  #geenide asukoht
  RNA_peak_plot = plotTranscripts(RNA_peak_annot$peak_list, RNA_peak_annot$peak_list, RNA_peak_annot$peak_annot, rescale_introns = FALSE, 
                                  region_coords = region_coords, connect_exons = FALSE, transcript_label = FALSE)
}
#/////////////////

#Construct metadata df for wiggleplotr
RNA_track_data = wiggleplotrGenotypeColourGroup(RNA_meta_df, gene_id, vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "RNA-seq")

RNA_coverage = plotCoverage(exons = CTCF_peak_annot$peak_list, cdss = CTCF_peak_annot$peak_list, track_data = RNA_track_data, rescale_introns = FALSE, 
                            transcript_annotations = CTCF_peak_annot$peak_annot, fill_palette = getGenotypePalette(),
                            connect_exons = FALSE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                            region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "line")

rm(RNA_counts)

#----------------------------------------------------------------------------------------------#
#'
#'ATAC seq andmete sisse lugemine
#'

#proovid
ATAC_sample_metadata = read.table(paste0(path,'QC_measures/run_sample_accession_PhaseIII.txt'))
colnames(ATAC_sample_metadata) = c("sample_id", "genotype_id")
ATAC_sample_metadata["condition_name"] = rep("naive", dim(ATAC_sample_metadata)[1])


#genotüübid
vcf_file = readRDS(paste0(path, 'genotypes/Kumasaka_100_samples.merged.rds'))
#ATAC_variant_information = importVariantInformation(paste0(path,"genotypes/Kumasaka_variantInformation.txt.gz"))


#Visualise caQTLs from the same loci
# ATAC seq lugemite arv used for calculating library size
ATAC_counts = read.table(paste0(path, "counts/counts_matrix.bed"), header = TRUE)



# atac metadata
ATAC_peak_metadata = ATAC_counts[,c("pid","start", "end","Chr", "strand")]
colnames(ATAC_peak_metadata)[4] = "chr"
ATAC_meta_df = wiggleplotrConstructMetadata(ATAC_counts[,7:dim(ATAC_counts)[2]], ATAC_sample_metadata, paste0(path, "bigwig"), condition_name_levels = c("naive"))


#Fetch all peaks in the region and make peak annot plot
ATAC_peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 7, ATAC_peak_metadata)

#////////////////
if (FALSE){
  ATAC_peak_plot = plotTranscripts(ATAC_peak_annot$peak_list, ATAC_peak_annot$peak_list, ATAC_peak_annot$peak_annot, rescale_introns = FALSE, 
                                   region_coords = region_coords, connect_exons = FALSE, transcript_label = FALSE)
}
#////////////////

#Construct metadata df for wiggleplotr
ATAC_track_data = wiggleplotrGenotypeColourGroup(ATAC_meta_df, rs_id, vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "ATAC-seq")

ATAC_coverage = plotCoverage(exons = ATAC_peak_annot$peak_list, cdss = ATAC_peak_annot$peak_list, track_data = ATAC_track_data, rescale_introns = FALSE, 
                             transcript_annotations = ATAC_peak_annot$peak_annot, fill_palette = getGenotypePalette(),
                             connect_exons = FALSE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "line")


rm(vcf_file,ATAC_counts)



#----------------------------------------------------------------------------------------------#
if (FALSE){
  #koik kuhmud, millel leidub significant QTL
  ca = read.table(paste(path, "tabix/cqn_permutations_100000_significant.sorted.txt.gz", sep=""))
  peak_ranges = GRanges(seqnames = ca$V2, strand = c("*"), ranges = IRanges(start = ca$V3, end = ca$V4, names = ca$V1))
  
  # kohad, kus kromatiin on avatud, leidub avatud kromatiiniga seotud statistiliselt oluline caQTL, sellega R2 seotud splaissimine ja CTCF seondumine
  interesting_peaks = subsetByOverlaps(peak_ranges, ctcf_contained_ranges)
  
  for (i in 1:length(interesting_peaks)){
    start = start(ranges(interesting_peaks))[i]
    end = end(ranges(interesting_peaks))[i]
    ca_7 = ca[which(ca$V10>start-3000 & ca$V10<end+3000),]
    if (dim(ca_7)[1]!=0){print(ca_7)}
  }
  
  
  #Make manhattan plot of the caQTL p-values
  ca_7 = ca[which(ca$V2==7 & ca$V10>region_coords[1] & ca$V10<region_coords[2]),]
  # Manhattan plot koigist caQTLdest mis asuvad 7 kromosoomil sellises vahemikus kus leidub moni huvipakkuv ctcf qtl
  library(ggplot2)
  peak_manhattan = ggplot(data = ca_7, aes(x = ca_7$V10, y = -log(ca_7$V19, 10))) +
    geom_jitter() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +  xlim(region_coords[1],region_coords[2])
}


#----------------------------------------------------------------------------------------------#


# formatting
ATAC_coverage$coverage_plot = ATAC_coverage$coverage_plot + labs(x='CQN')


RNA_coverage$coverage_plot = RNA_coverage$coverage_plot + labs(x='CQN')
ATAC_coverage$coverage_plot$theme$axis.ticks.x = element_blank()


#Make a joint plot
joint_plot = cowplot::plot_grid(ATAC_coverage$coverage_plot + labs(x='CQN'), 
                                RNA_coverage$coverage_plot + labs(x='CQN'),
                                ATAC_coverage$tx_structure + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()), 
                                RNA_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(3,3,2,2))


library(ggplot2)
ggsave(paste0("/gpfs/hpchome/evelin95/plots/",rs_id,".pdf"), plot = joint_plot, width = 3.5, height = 6)



