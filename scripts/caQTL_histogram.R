
path = "/gpfs/hpchome/a72094/rocket/projects/chromatin_to_splicing/QTL/"
caQTL = read.table(paste0(path,"FDR_corrected/cqn_permutations_100000_significant.txt"))


within_peak <- function(peak){
  if (peak[10] < peak[4] & peak[3] < peak[10]) return(as.factor(1))
  else return(as.factor(0))
}

distance <- function(peak){
  return(as.numeric(peak[3]) - as.numeric(peak[10]) + (as.numeric(peak[4]) - as.numeric(peak[3]))/2 )
}

caQTL$within_peak = apply(caQTL, 1, within_peak)
caQTL$dist = apply(caQTL, 1, distance)

dim(caQTL)
dim(caQTL[which(abs(caQTL$dist) < 1000),])
dim(caQTL[which(caQTL$within_peak==1), ])


library(ggplot2)



sub = as.data.frame(caQTL[which(abs(caQTL$dist) < 1000),c("within_peak", "dist")])
ggplot(sub, aes(x = dist, colour = within_peak)) + 
  geom_histogram(binwidth = 10, fill = "white") +
  labs(x = "caQTL kaugus kühmu keskpunktist", y = "Sagedus") +
  guides(color = "legend") +
  scale_color_manual(labels = c("väljaspool avatud piirkonda", "avatud piirkonna sees"), values = c("royalblue2", "tomato2")) +
  theme(legend.position = c(0.78, 0.92), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "white",size = 2, linetype = "solid", color = "grey80"),   
        panel.grid.major = element_line(size = 0.5, linetype = 'solid'), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "grey80")) 

ggsave("/gpfs/hpchome/evelin95/plots/ca_histogramm.pdf",  width = 6, height = 5)
