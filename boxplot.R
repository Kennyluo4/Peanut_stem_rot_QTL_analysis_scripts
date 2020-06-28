# plot Boxplot of traits per QTL
dt = read.csv('QTLpeakdata.csv')
relation = read.csv('peak_marker.csv')
dt[,2:30] = lapply(dt[,2:30],as.numeric)
for (i in 1:104){  #i is the number of QTL in "peak_marker.csv", loop over every QTL/in each row
  #browser()
  markername = as.character(relation[i,2])
  envname = as.character(relation[i,1])
  QTLname = as.character(relation[i,3])
  PICname = paste(QTLname,envname,'tiff',sep = '.')
  modelname = lm(dt[,envname]~dt[,markername], data = dt)    #(environment phenotype~ marker genotype)
  result = anova(modelname)
  stat_p = ifelse(result$`Pr(>F)`<= 0.001, "<0.001", round(result$`Pr(>F)`,digits = 3))
  #stat_p = round(result$`Pr(>F)`,digits = 3)
  tiff(filename = PICname,width = 1200, height = 1000, res = 200)
  boxplot(dt[,envname]~dt[,markername], main = paste(QTLname,"\n",envname), ylab = 'phenotype')
  #legend("topright", bty="n", legend=paste('p-value:', stat_p[1]))
  title(sub = paste('p-value:', stat_p[1]))
  dev.off()
  }

