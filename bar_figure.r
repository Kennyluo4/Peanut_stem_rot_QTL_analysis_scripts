# plot bar figure of genes' functions withn QTLs
# input file: 1st column as QTL, 2nd column genes' functions with QTL regions
library(RColorBrewer)
library(ggplot2)
mydt = read.csv('test.csv')

tiff("Fig 2.tiff", width = 4500, height = 2600, res = 500)
  ggplot(mydt) +
    aes(x = qtl,fill=class) +
    #scale_fill_brewer(palette="Set1") +
    #scale_color_gradientn(colours = rainbow(12)) +
    geom_bar() +
    labs(x = "", y = "Frequency") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=7))
dev.off()
