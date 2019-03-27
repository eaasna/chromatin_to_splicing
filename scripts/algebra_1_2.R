

lahendid <- data.frame(
  r = c(1, 1, 1, 1, 2.45, 2.45, 2.45, 2.45),
  phi = c(0, 90, 180, 270, 67.5, 157.5, 247.5, 347.5)
)

x = round(lahendid$r * cos(lahendid$phi/180*pi),2)
y = round(lahendid$r * sin(lahendid$phi/180*pi),2)
koordinaadid = paste(x,y,sep=",")

punktid = data.frame(x, y)

library(ggplot2)
  ggplot(punktid,aes(x,y))+
  geom_point(aes()) + 
  labs(x="reaaltelg", y="imaginaartelg") +
  geom_label(aes(x+.1,y+.1,label=koordinaadid))
