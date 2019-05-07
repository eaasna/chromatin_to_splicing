library(readr)
library(Rsamtools)
library(ggplot2)

path = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/"

# koik atac kuhmud millel on qtl
ca = read.table(paste(path, "tabix/cqn_permutations_100000_significant.sorted.txt.gz", sep=""))
ca_peaks = GRanges(seqnames = ca$V9, strand = c("*"), ranges = IRanges(start = ca$V10, end = ca$V11, names = ca$V8))


#koik ctcf seondumiskohad
ctcf = read.table("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/Blood_ATAC/processed/CTCF/qtltools/output/cqn/CTCF.permuted.txt.gz")
bindingsites <- GRanges(seqnames = ctcf$V2, strand = c("*"), ranges = IRanges(start = ctcf$V3, end = ctcf$V4, names = ctcf$V1))


# kuhmud mis paiknevad ctcf seondumiskohas
ca_peaks_with_binding_site = subsetByOverlaps(ca_peaks, bindingsites)
ca_within_peaks = ca[which(ca$V8 %in% names(ca_peaks_with_binding_site)), ]

# QTLd mis seotud kromatiini avatuse ja CTCF seondumisega
# V4 on ctcf QTLd
r2 = read.table(paste(path, "results/ctcf/cqn_rsq08.txt",sep=""))
# kontroll
ctcfQTL = read.table(paste(path,"CTCF/CTCF.permuted.significant.sorted.txt.gz", sep = ""))
ctcfQTL[which(ctcfQTL$V8 %in% r2$V4),]

# kuhmud mille qtl seotud ctcf qtl'ga
ca_r2 = ca[which(ca$V8 %in% r2$V1), ]

# vaja teha 1 ggplot, millel kujutatud erinevate varvidega 3 pvalue histogrammi
# pvaartused failidest ca, ca_within_peaks, ca_r2

#stacking 3 dataframes
ca_values<-ca[,"V19"]
ca_withinpeak_values<-ca_within_peaks[,"V19"]
ca_r2_values<-ca_r2[,"V19"]
all<-data.frame(dataset=c(rep('all',length(ca_values)),rep('withinpeak',length(ca_withinpeak_values)), rep('r2',length(ca_r2_values))),value=c(ca_values,ca_withinpeak_values, ca_r2_values))


ggplot(all,aes(x=value,fill=dataset))+
  geom_histogram(aes(y=0.001*..density..),
                 alpha=0.3,position='identity', binwidth = 0.001) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')



ggplot(all[which(all$dataset=='all'),],aes(x=value,fill=dataset))+
  geom_density(aes(y=0.01*..density..), alpha=0.3) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')

ggplot(all[which(all$dataset=='withinpeak'),],aes(x=value,fill=dataset))+
  geom_density(aes(y=0.01*..density..), alpha=0.3) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')

ggplot(all[which(all$dataset=='r2'),],aes(x=value,fill=dataset))+
  geom_density(aes(y=0.01*..density..), alpha=0.3) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')


ggplot(all,aes(x=value,fill=dataset))+
  geom_density(aes(y=0.01*..density..), alpha=0.3) +
  facet_wrap(~dataset,nrow=3) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')


ggplot(all,aes(x=value,fill=dataset))+
  geom_histogram(aes(y=0.001*..density..), binwidth = 0.001)+
  geom_density(aes(y=0.001*..density..), alpha=0.3) +
  facet_wrap(~dataset,nrow=3) +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')


ggplot(all[which(all$dataset=='all'),], aes(x=value))+
  geom_histogram(aes(y=0.01*..density..), binwidth = 0.01)+
  geom_density() +
  labs(x='p-väärtus', y='tihedus', title='Tihedus ehk probability density')
