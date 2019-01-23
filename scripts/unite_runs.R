data = read.table("sample_accession_run_accession.txt", header=TRUE)
path = "/gpfs/hpchome/a72094/rocket/datasets/Kumasaka_2017_ATAC/PRJEB9977/"
outpath = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/merged_fastq/"
head(data)
library("dplyr")
samples = unique(data$secondary_sample_accession)

sample_runs = dplyr::select(dplyr::filter(data, data$secondary_sample_accession == samples[1]), run_accession)

#loopin läbi kõik samples
#leian iga sampliga kaasnevad runid 
#catin iga sample kohta kaks faili
#*sample_name*_1.fastq.gz
#*sample_name*_2.fastq.gz
#loopin läbi kõik sample_runs ja koostan kaks faili


#Looking at all samples
for (samplename in samples) {
  sample_runs = dplyr::select(dplyr::filter(data, data$secondary_sample_accession == samplename), run_accession)
  cat(paste("cat ", sep = ""), file = "cat_1.txt", append = TRUE)
  cat(paste("cat ", sep = ""), file = "cat_2.txt", append = TRUE)
  for (runname in sample_runs) {
    cat(paste("/gpfs/hpchome/a72094/rocket/datasets/Kumasaka_2017_ATAC/PRJEB9977/",runname,"/",runname, "_1.fastq.gz", sep = ""), file = "cat_1.txt", append = TRUE)
    cat(paste("/gpfs/hpchome/a72094/rocket/datasets/Kumasaka_2017_ATAC/PRJEB9977/",runname,"/",runname, "_2.fastq.gz", sep = ""), file = "cat_2.txt", append = TRUE)
  }
  cat(paste(" > /gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/merged_fastq/", samplename, "_1.fastq.gz \n", sep = ""), file = "cat_1.txt", append = TRUE)
  cat(paste(" > /gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/merged_fastq/", samplename, "_2.fastq.gz \n", sep = ""), file = "cat_2.txt", append = TRUE)
}
