map<-read.csv(file = "mapinfomation.csv",header = TRUE,as.is=TRUE)
QTL<-read.csv(file = "QTLinfomation.csv",header = TRUE,as.is=TRUE)
library(circlize)
tiff("Fig 2.tiff", width = 3600,height = 3580,res = 600)
par(mar = c(0.01, 0.01, 0.01, 0.01), lwd = 0.3, cex = 0.7)
#circos.par("track.height" = 0.1)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
#dinitialize the map layout, assign factors and locations
circos.initialize(factor=map$LinkAgegroup,x = map$LociPosition)
#create a rigion/circle for the plot, specify the height of the regions, x-axels as linkage group coordinate
circos.trackPlotRegion(factors =map$LinkAgegroup,track.height = 0.14,ylim = c(-1.2,1.2),
panel.fun = function(x, y) {
circos.axis(major.at = seq(from = 0,to = 300,by = 20),labels = paste(seq(0,300,by = 20),sep = ),
lwd = 0.7, major.tick.percentage=0.4,labels.away.percentage=0.05,
minor.ticks = 10,labels.cex = 0.4)
})

#plot the markers position on the region just created
for(i in 1:length(map[,1])){
circos.lines(x = c(map[i,3],map[i,3]),y = c(-1.52,1.51),sector.index =map[i,2],track.index = 1,col = "light blue",type = "l",lwd = 0.6,lty = 1)
}

#create another rigion/circole for ploting the QTLs LOD scores
circos.trackPlotRegion(factors =map$LinkAgegroup,track.height = 0.14, ylim = c(0,12.8))
circos.trackLines(factor=QTL$Chromosome,x = QTL$Peak,y = QTL$lodscore,col = "orange",lwd = 1.2 ,straight = TRUE,type = "h",area = TRUE)

#plot the QTLs 
for(m in 1:length(QTL[,1])){
  circos.lines(x = c(QTL[m,7],QTL[m,8]),y = c(-6, -6),sector.index = QTL[m,2],track.index = 1,lwd = 1.4, col = "red")
}
for(m in 1:length(QTL[,1])){
  circos.lines(x = c(QTL[m,6],QTL[m,6]),y = c(-5.8,-6), sector.index = QTL[m,2],track.index = 1, lwd = 1, col = 'red')
}

#epistatic effect qSR.A01-1 - qSR.B07_1-1
#circos.link(, rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A05/B05',c(154.011,162.364),'A05/B05',c(168.37,171.017), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('B02',c(192.484,202.012),'B02',c(202.012,210.84), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A08',c(23.782,43.266),'A10/B04',c(87.349,92.687), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A03/B03',c(28.346,30.396),'A05/B05',c(147.419,154.011), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('B01',c(107.622,123.462),'B02',c(223.824,225.04), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A01',c(8.52900000000001,25.1),'B05',c(21.876,26.301), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A05/B05',c(89.823,91.177),'B06_1',c(126.291,130.126), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A02',c(52.638,55.769),'B09',c(19.852,22.315), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('B03',c(6.145,11.786),'B09',c(44.603,51.332), rou1 = 0.645,rou2 = 0.645, col = "light green")
circos.link('A01',c(51.09,51.91),'B07_1',c(10.6,10.6), rou1 = 0.645,rou2 = 0.645, col = "light green")

#on circle region1, add the chromosome name
circos.trackPlotRegion(factors = map$LinkAgegroup,track.index = 1,
                       panel.fun = function(x, y) {
                         sector.index = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), -0.1, sector.index,cex = 0.8,font = 2,facing = 'clockwise')
                       })


dev.off()

#circos.initializeWithIdeogram(cytoband = system.file(package = "circlize",
#                                                     "extdata", "cytoBand.txt"), species = NULL, sort.chr = TRUE,
 #                             chromosome.index = NULL, major.by = NULL,
 #                             plotType = c("ideogram", "axis", "labels"),
  #                            track.height = NULL, ideogram.height = convert_height(2, "mm"),
   #                           )
#circos.track(ylim = c(0, 1))
#circos.genomicIdeogram()
#read.cytoband(cytoband = system.file(package = "circlize",
#                                     "extdata", "cytoBand.txt"), species = NULL, chromosome.index = NULL, sort.chr = TRUE)

