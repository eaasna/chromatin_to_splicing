#Secondary sample accession data currently used for naming bam files
#siin on 100 proovi, 91 neist on Phase III. Teised jätan kõrvale 
samples_list = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QC_measures/sample_accession.txt")

sink("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QC_measures/run_sample_accession.txt")
for (i in 1:length(samples_list[,1])){
  sample = as.character(samples_list[i,1])
  sample = "ERS798625"
  bamstat = read.table(paste("/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/bamstat/",sample, ".bamstat.txt", sep =""), header = TRUE)

  full = dplyr::transmute(bamstat, SampleID = bamstat$SampleID, hom_consistency=bamstat$n_hom_consistent/bamstat$n_hom_covered, het_consistency=bamstat$n_het_consistent/bamstat$n_het_covered)
  full = dplyr::mutate(full, total_consistency=full$hom_consistency+full$het_consistency)

  #sampleID used in vcf file to reference genotypes
  #max_rows=dplyr::filter(full, full$total_consistency >= 0.95*max(full$total_consistency))
  max_rows=dplyr::filter(full, full$total_consistency >= 0.95*2)
  sampleID = as.character(max_rows[1,1])
  #print(max_rows)
  
  cat(paste(sample, sampleID, "\n", sep="\t"))
  
  #Plotting
  library(ggplot2)
  hom = full[which(full[,"SampleID"]==sampleID),]$hom_consistency
  het = full[which(full[,"SampleID"]==sampleID),]$het_consistency
  hom_other = full[which(full[,"SampleID"]!=sampleID),]$hom_consistency
  het_other = full[which(full[,"SampleID"]!=sampleID),]$het_consistency
  plot = ggplot(data = full[which(full[,"SampleID"]!=sampleID),]) + 
    geom_point(mapping = aes(het_other, hom_other), color = "tomato2") + 
    geom_point(mapping = aes(het,hom), color = "green3") +
    labs(x="Kokkulangevus heterosügootides", y="Kokkulangevus homosügootides") + 
    coord_fixed(ratio = 1) +
    annotate(geom = "text", x = het - 0.03, y = hom - 0.03, label = sampleID) +
    annotate(geom = "text", x = mean(het_other) - 0.03, y = mean(hom_other) - 0.03, label = "Sobimatud genotüübid") +
    expand_limits(x = 0, y = 0) +
    theme(panel.background = element_rect(fill = "white", colour = "white",size = 2, linetype = "solid", color = "grey80"),   
          panel.grid.major = element_line(size = 0.5, linetype = 'solid'), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "grey80"))
  
  print(plot)
  
  ggsave("/gpfs/hpchome/evelin95/plots/sample_genotype_match.pdf",  width = 5, height = 5)

}

sink()






  