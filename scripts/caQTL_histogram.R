caQTL = read.table("\\Users\\evelin95\\Documents\\chromatin_to_splicing\\cqn_permutations_100000_significant.sorted.txt.gz")


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

sub = as.data.frame(caQTL[which(abs(caQTL$dist) < 10000),c("within_peak", "dist")])
ggplot(sub, aes(x = dist, colour = within_peak)) + 
  geom_histogram(binwidth = 100, fill = "white" ) +
  labs(x = "caQTL kaugus kühmu keskpunktist", y = "Sagedus") +
  guides(color = "legend") +
  scale_color_discrete(labels = c("väljaspool avatud piirkonda", "avatud piirkonna sees")) +
  theme(legend.position = c(0.78, 0.92), legend.title = element_blank())

ggsave("\\Users\\evelin95\\Documents\\chromatin_to_splicing\\plots\\caQTL_histogram.png", width=5, height = 5)
