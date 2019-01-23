library(readr)
library(Rsamtools)
library(dplyr)
load_all("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/seqUtils/")

add_to_df_with_info <- function(rida) {
  e_id = e_names[rida['expression']]  #snp_id from genotypes file
  ca_id = ca_names[rida['chromatin_accessibility']] 
  gene_line = eQTL[[1]][which(eQTL[[1]]$genotype_snp_id == e_id),c(1, 2, 3, 4, 6, 8, 10, 16, 19)]
  peak_line = caQTL[[1]][which(caQTL[[1]]$genotype_snp_id == ca_id),c(1, 2, 3, 4, 6, 8, 10, 12, 13)]
  colnames(gene_line) <- c("gene_id", "gene_chr", "gene_start", "gene_end", "gene_n_snps", "gene_snp_id", "gene_snp_location", "gene_p_nominal", "gene_beta")
  colnames(peak_line) <- c("peak_id", "peak_chr", "peak_start", "peak_end", "peak_n_snps", "peak_snp_id", "peak_snp_location", "peak_p_nominal", "peak_beta")
  
  rsq = rsq_matrix[rida['chromatin_accessibility'], rida['expression']]
  gene_line = dplyr::mutate(gene_line, rsquared = rsq)
  
  if (nrow(peak_line)==nrow(gene_line)){        #!! peab tegelema mitmerealistega
    pair = dplyr::bind_cols(peak_line, gene_line)
  } else {
    pair = dplyr::bind_cols(peak_line[1,], gene_line[1,])
  }
  
  write.table(pair, file = outfile, col.names = F, quote = F, row.names = F, append = T)
}

add_to_df_coverage <- function(rida) {
  e_id = e_names[rida['expression']]  #snp_id from genotypes file
  genes = eQTL[[1]][which(eQTL[[1]]$genotype_snp_id == e_id),c(1, 2, 3, 4, 6, 8, 10, 16, 19)]
  colnames(genes) <- c("gene_id", "gene_chr", "gene_start", "gene_end", "gene_n_snps", "gene_snp_id", "gene_snp_location", "gene_p_nominal", "gene_beta")
  write.table(genes, file = outfile, col.names = F, quote = F, row.names = F, append = T)
}

add_to_df <- function(rida) {
  e_id = e_names[rida['expression']]  #snp_id from genotypes file
  ca_id = ca_names[rida['chromatin_accessibility']] 
  if (e_id != ca_id) {
    gene_line = eQTL[[1]][which(eQTL[[1]]$genotype_snp_id == e_id),c(8, 9, 10)]
    peak_line = caQTL[[1]][which(caQTL[[1]]$genotype_snp_id == ca_id),c(8, 9, 10)]
    colnames(gene_line) <- c("gene_snp_id", "gene_snp_chr", "gene_snp_location")
    colnames(peak_line) <- c("peak_snp_id", "peak_snp_chr", "peak_snp_location")
    
    rsq = rsq_matrix[rida['chromatin_accessibility'], rida['expression']]
    gene_line = dplyr::mutate(gene_line, rsquared = rsq)
    
    pair = dplyr::bind_cols(peak_line[1,], gene_line[1,])
    
    write.table(pair, file = outfile, col.names = F, quote = F, row.names = F, append = T)
    
  }
}

#sub = "upstream"
#sub = "contained"
#sub = "downstream"
sub = "featureCounts"

# overap between ATAC peaks and RNA QTLs that is covered peaks
#outfile = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/rsquared/cqn_txrevise_", sub, "_coverage.txt", sep ="")
#eqtl_file = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/txrevise.significant.", sub, ".sorted.txt.gz", sep="")


# overlap between covered ATAC peaks and CTCF QTLs 
#outfile = paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_", sub, "_rsq.txt", sep="")
#eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.sorted.txt.gz"

# overlap between all ATAC peaks and CTCF QTLs 
outfile = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/results/ctcf/cqn_rsq.txt"
eqtl_file = "/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/CTCF/CTCF.permuted.significant.sorted.txt.gz"


#-----------

genotype = readRDS("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/genotypes/Kumasaka_100_samples.merged.rds")

rsq <- function (x, y) cor(x, y) ^ 2

chr_lengths = read.table("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/chr_lengths.txt", header = F)
colnames(chr_lengths) <- c("chr", "length")

chrs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
for (chr in chrs){
  print(chr)
  start = 1
  len = chr_lengths[chr,2]
  pairs = data.frame(peak_id=as.character(), stringsAsFactors=FALSE) #tühi dataframe, millesse kogun paarid
  
#2e+06 nii suures aknas otsin variante
  while (start < len - 2e+06) {
    gr <- GRanges(seqnames = chr, strand = c("*"),
                  ranges = IRanges(start = start, width=2e+06))
    start = start + 1e+06
    

    eQTL = scanTabixDataFrame(eqtl_file, gr, col_names = FALSE)

    # finding covered peaks
    caQTL = scanTabixDataFrame("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_100000_significant.sorted.txt.gz", gr, col_names = FALSE)
    
    # finding covered peaks that are enriched with CTCF
    #caQTL = scanTabixDataFrame(paste("/gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_", sub, "_enriched.sorted.txt.gz", sep=""), gr, col_names = FALSE)
    
    if (lengths(eQTL[])[1]!=0 & lengths(eQTL[])[1]!=1 & lengths(caQTL[])[1]!=0 & lengths(caQTL[])[1]!=1){
      
      # !! ainult 3 tüübiga QTLtools väljundi jaoks <- 21 tulpa
      #eQTL[[1]]=eQTL[[1]][,c(-6,-7)]
      colnames(eQTL[[1]]) <- c("phenotype_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "n_snps", "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "freedom", "dummy", "first", "second", "p_nominal", "slope", "p_empirical", "p_beta")
      
      colnames(caQTL[[1]]) <- c("phenotype_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "n_snps", "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "freedom", "dummy", "first", "second", "p_nominal", "slope", "p_empirical", "p_beta")
      #eemalda caQTL ja eQTL andmetest sellised lead SNP'd mille asukohta pole genotüübitud
      
      snppos = unlist(genotype$snpspos[which(genotype$snpspos$chr == chr), 'pos'], use.names = F)
      eQTL[[1]] = dplyr::filter(eQTL[[1]], eQTL[[1]]$snp_start %in% snppos) #ainult ühe kromosoomi andmed
      caQTL[[1]] = dplyr::filter(caQTL[[1]], caQTL[[1]]$snp_start %in% snppos) #ainult ühe kromosoomi andmed
      
      #QTL ja genotypes andmete snp_id'd ei kattu
      #küsida QTL snp start järgi 
      e_snps_pos = eQTL[[1]][which(eQTL[[1]]$pheno_chr == chr), 'snp_start']
      e_snps = unlist(dplyr::filter(genotype$snpspos, genotype$snpspos$pos %in% unlist(e_snps_pos))[,'snpid'], use.names = F)
      e_genotypes = subset(genotype$genotypes, rownames(genotype$genotypes) %in% e_snps)
      
      #lisada eQTL tabelisse genotype faili vastavad snp_id'd <- see on enamasti sama nagu QTL failide snp_id samal positsioonil, aga mõnel juhul erineb
      geno_osa = dplyr::filter(genotype$snpspos, genotype$snpspos$pos %in% unlist(e_snps_pos))
      geno_osa = geno_osa[which(geno_osa$chr == chr), c('pos', 'snpid')]
      eQTL[[1]] = dplyr::left_join(eQTL[[1]], geno_osa, by = c('snp_start'='pos'))
      colnames(eQTL[[1]])[20] <- "genotype_snp_id"
      
      ca_snps_pos = caQTL[[1]][which(caQTL[[1]]$pheno_chr == chr), 'snp_start']
      ca_snps = unlist(dplyr::filter(genotype$snpspos, genotype$snpspos$pos %in% unlist(ca_snps_pos))[,'snpid'], use.names = F)
      ca_genotypes = subset(genotype$genotypes, rownames(genotype$genotypes) %in% ca_snps)
      
      #lisada caQTL tabelisse genotype faili vastavad snp_id'd
      geno_osa = dplyr::filter(genotype$snpspos, genotype$snpspos$pos %in% unlist(ca_snps_pos))
      geno_osa = geno_osa[which(geno_osa$chr == as.character(chr)), c('pos', 'snpid')]
      caQTL[[1]] = dplyr::left_join(caQTL[[1]], geno_osa, by = c('snp_start'='pos'))
      colnames(caQTL[[1]])[20] <- "genotype_snp_id"
      
      #kas selles piirkonnas leidub SNP'e
      if (length(e_genotypes)!=0 & length(ca_genotypes)!=0){
        # mul on kaks maatriksit ja tahan leida iga nende kahe rea kombinatsiooni cor ruudu
        # rows from caQTL
        # columns from eQTL
        rsq_matrix = apply(e_genotypes, 1, function(rida) apply(ca_genotypes, 1, rsq, y = rida)) #!! kui ca või e snp'e on ainult 1 tekib vektor mitte maatriks
        
        
        #kas ca_genotype ja e_genotype hulgas on rohkem kui 1 snp id
        
        
        if (is.matrix(rsq_matrix)){
          if (dim(rsq_matrix)[2]!=1){
            rsq_matrix = rsq_matrix[complete.cases(rsq_matrix), ]
          }
        } else {
            rsq_matrix = matrix(c(rsq(ca_genotypes[1,], e_genotypes[1,])))
            #seda saaks ilusamaks teha apply kasutades
            if (dim(e_genotypes)[1]>1){
              for (i in 2:dim(e_genotypes)[1]){
                rsq_matrix = cbind(rsq_matrix, c(rsq(ca_genotypes[1,], e_genotypes[i,])))
              }
            }
            rownames(rsq_matrix) <- row.names(ca_genotypes)
            colnames(rsq_matrix) <- row.names(e_genotypes) #!!! viga
        }
        #kuidas leida indeks, kus suurem kui 0.8
        significant = as.matrix(which(rsq_matrix > 0.8, arr.ind = T))
        
        #kas leidub mõni paar, mille rsquared>0.8
        if (dim(significant)[1]!=0){
          row.names(significant) <- NULL
          colnames(significant) <- c("chromatin_accessibility", "expression")
          
          ca_names = rownames(rsq_matrix)
          e_names = colnames(rsq_matrix)
          
          apply(significant, 1, add_to_df_coverage)
        }
      }
      
    } #// if 
    
  } #// while
  
  
  #!!!!!
  pairs = read.table(outfile)
  write.table(unique(pairs), file = outfile, col.names = F, quote = F, row.names = F, append = F)
}

#piirjuhtude handlimine saaks kindlasti paremini lahendada, töötab, aga on kole vaadata
#kirjutada eraldi funktsioon, mida applida ja kasutada selle sees rsq funktsiooni? 

#leia kõik geenid, millel on vähemalt üks kühm
#kõik geenid mille põhjuslik variant on jagatud mingi kühmuga, kühmud jätan kirjutamata?