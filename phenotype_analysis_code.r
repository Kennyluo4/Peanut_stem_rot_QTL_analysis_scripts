# draw the phenotype distribution
library(ggplot2)
library(gridExtra)
dt = read.csv('combined_average.csv', row.names = 1)
dt2 = read.csv('combined_average_class.csv')
#########version1 individual experiment plots
pdf('Data_distribution.pdf')
par(mfrow=c(3,4))
for(i in 1:ncol(dt)){
  hist(dt[,i],prob=T,col="light blue", main = colnames(dt)[i],xlab = "", breaks = 15)
  xfit<-seq(min(dt[,i], na.rm = TRUE),max(dt[,i], na.rm = TRUE),length=40)
  yfit<-dnorm(xfit,mean(dt[,i], na.rm = TRUE),sd(dt[,i],na.rm = TRUE))
  lines(xfit,yfit,col="red",lwd=3)
}
dev.off()

########version 2 individual experiment plot using ggplot2
tiff('Data_distribution_test.tiff',width = 4300, height = 3200, res= 450)
p <- list()
for (i in colnames(dt)) {
  p[[i]] = ggplot(data = dt, aes_string(x = dt[,i])) +
    geom_histogram(bins = 18, alpha = 0.8, color="darkblue", fill="lightblue", size =0.1) + 
#    geom_text(data=dt, aes(label="Tifrunner", x=dt["Tifrunner",i], y= 2+2), size=4) +
#    geom_segment(data=dt, 
#                 aes(x=dt["Tifrunner",i], xend=dt["Tifrunner",i], y=2 + 1.5, yend=2 + 0.25),
#                 arrow=arrow(length=unit(2, "mm")), colour = "blue") +
#    geom_text(data=dt, aes(label="NC3033", x=dt["NC 3033",i], y= 2+2), size=4) +
#    geom_segment(data=dt, 
#                 aes(x=dt["NC 3033",i], xend=dt["NC 3033",i], y=2 + 1.5, yend=2 + 0.25),
#                 arrow=arrow(length=unit(2, "mm")), colour = "green") +
    labs(x = i) + ylab("")
}
do.call(grid.arrange,p)
dev.off()

#test
ggplot(data = dt, aes(x = dt[,"sc_13"])) +
  geom_histogram(bins = 18, alpha = 0.8, color="darkblue", fill="lightblue", size =0.1) + 
  geom_text(data=dt, aes(label="Tifrunner", x=dt["Tifrunner",i], y= 2+2), size=4) +
  geom_segment(data=dt, 
               aes(x=dt["Tifrunner",i], xend=dt["Tifrunner",i], y=2 + 1.5, yend=2 + 0.25),
               arrow=arrow(length=unit(2, "mm")), colour = "blue") +
  geom_text(data=dt, aes(label="NC3033", x=dt["NC 3033",i], y= 2+2), size=4) +
  geom_segment(data=dt, 
               aes(x=dt["NC 3033",i], xend=dt["NC 3033",i], y=2 + 1.5, yend=2 + 0.25),
               arrow=arrow(length=unit(2, "mm")), colour = "green") +
  labs(x = i) + ylab("")





#######plot distributions of multiple experiment 

library(ggplot2)
library(gridExtra)
#pdf('Data_distribution_1.pdf',width = 11, height = 8)
tiff('Data_distribution_1.tiff',width = 4700, height = 3500, res= 600)
plot1 = ggplot(data = subset(dt2,Trait =="score"),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Disease score") +
  labs(tag = "a")
plot2 = ggplot(data = subset(dt2,Trait =="percentage"),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Percentage of resistant plants") +
  labs(tag = "b")
plot3 = ggplot(data = subset(dt2,(Trait =="broadcast"&Year!= "13")),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Hits/foot") +
  labs(tag = "c")
plot4 = ggplot(data = subset(dt2,(Trait =="BF")),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Disease dispersal index") +
  labs(tag = "d")
plot5 = ggplot(data = subset(dt2,(Trait =="WF"&Year!= "13")),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Disease dispersal index") +
  labs(tag = "e")
#incorporate BF13 into plot4
#plot6 = ggplot(data = subset(dt2,Experiment =="BF13"), aes(x = Phenotype,fill = Experiment)) +
#  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Disease incidence")
plot7 = ggplot(data = subset(dt2,(Experiment ==c("BA13", "BU13"))),aes(x = Phenotype,fill = Experiment)) +
  geom_histogram(bins = 20, alpha = 0.68, position = "dodge") + labs(x = "Disease severity") +
  labs(tag = "f")
grid.arrange(plot1,plot2,plot3,plot4,plot5,plot7)
dev.off()


##normality test
normal.test = sapply(dt,function(i) shapiro.test(i)$p.value)
# QQ plot
pdf('QQ-plot.pdf')
par(mfrow=c(3,4))
for(i in 1:ncol(dt)){
  qqnorm(dt[,i], cex = 0.5, main=colnames(dt)[i], xlab = '')
  qqline(dt[,i])}
dev.off()

install.packages("moments")
library(moments)
kurtosis(dt$BF15_1,na.rm = TRUE)
kurtosis(dt$WF15_1,na.rm = TRUE)
skewness(dt$BF15_1, na.rm = TRUE)
skewness(dt$WF15_1, na.rm = TRUE)

##anova
fl15 = read.csv('FL2015.csv')
kt_fl13_broad = kruskal.test(phenotype~rate_time, data = fl13_broad)
kt_fl14_broad_geno_rate = kruskal.test(phenotype~interaction(fl14_broad$genotype,fl14_broad$rating),data = fl14_broad)
kt_fl15point_geno_broad = kruskal.test(phenotype~rate_time, data = fl15_point)

fl14_15 = read.csv('FL2014_2015.csv')
fl14_15_point = subset(fl14_15,method !="broad")
str(fl14_15_point)
fl14_15_point$rep = as.factor(fl14_15_point$rep)
kt_flcombine_point_geno = kruskal.test(phenotype~genotype, data = fl14_15_point)

pairwise.wilcox.test(x = fl14_15_point$phenotype, g = fl14_15_point$rating, p.adjust.method = "BH")
library(FSA)
dunnTest(phenotype~rating, data = fl14_15_point, method = "bh")
ggboxplot(fl14_15,x = 'method', y = 'phenotype',color = 'method', xlab = 'Inoculum', ylab = 'Disease score')

fl15_1 = read.csv("FL2015_1.csv")
fl15_1$rep = as.factor(fl15_1$rep)
fl15_1_bl = subset(fl15_1,method=="blue")
fl15_1_wf = subset(fl15_1,method=="white")
ano_fl15_bf_1 = aov(phenotype~genotype,data = fl15_1_bl)
summary(ano_fl15_bf)
ano_fl15_wf_1 = aov(phenotype~genotype,data = fl15_1_wf)
summary(ano_fl15_wf_1)

###correlation
install.packages("agricolae")
library(agricolae)
cor.table = correlation(dt, method = "spearman")
write.csv(cor.table$correlation,"correlation.csv")
write.csv(cor.table$pvalue,"correlation_pvalue.csv")





###caulculating AUDPC######
library(plantbreeding)
reading.dates_14 = as.Date(c("2014-09-26","2014-10-20","2014-10-24"))
reading.dates_15 = as.Date(c("2015-09-29","2015-10-10"))
all = read.csv('AUDPC_data.csv')
b_14_rep1 = all[,1:4]
b_14_rep2 = all[,c(1,5:7)]
bf_14_rep1 = all[,c(1,8:10)]
bf_14_rep2 = all[,c(1,11:13)]
bf_14_rep3 = all[,c(1,14:16)]
bf_14_rep4 = all[,c(1,17:19)]
wf_14_rep1 = all[,c(1,20:22)]
wf_14_rep2 = all[,c(1,23:25)]
wf_14_rep3 = all[,c(1,26:28)]
wf_14_rep4 = all[,c(1,29:31)]
b_15_rep1 = all[,c(1,32:33)]
b_15_rep2 = all[,c(1,34:35)]
bf_15_rep1 = all[,c(1,36:37)]
bf_15_rep2 = all[,c(1,38:39)]
wf_15_rep1 = all[,c(1,40:41)]
wf_15_rep2 = all[,c(1,42:43)]

AUDPC_b_14_rep1 = AUDPC.cal(reading.dates_14,b_14_rep1)
AUDPC_b_14_rep2 = AUDPC.cal(reading.dates_14,b_14_rep2)
AUDPC_bf_14_rep1 = AUDPC.cal(reading.dates_14,bf_14_rep1)
AUDPC_bf_14_rep2 = AUDPC.cal(reading.dates_14,bf_14_rep2)
AUDPC_bf_14_rep3 = AUDPC.cal(reading.dates_14,bf_14_rep3)
AUDPC_bf_14_rep4 = AUDPC.cal(reading.dates_14,bf_14_rep4)
AUDPC_wf_14_rep1 = AUDPC.cal(reading.dates_14,wf_14_rep1)
AUDPC_wf_14_rep2 = AUDPC.cal(reading.dates_14,wf_14_rep2)
AUDPC_wf_14_rep3 = AUDPC.cal(reading.dates_14,wf_14_rep3)
AUDPC_wf_14_rep4 = AUDPC.cal(reading.dates_14,wf_14_rep4)
AUDPC_b_15_rep1 = AUDPC.cal(reading.dates_15,b_15_rep1)
AUDPC_b_15_rep2 = AUDPC.cal(reading.dates_15,b_15_rep2)
AUDPC_bf_15_rep1 = AUDPC.cal(reading.dates_15,bf_15_rep1)
AUDPC_bf_15_rep2 = AUDPC.cal(reading.dates_15,bf_15_rep2)
AUDPC_wf_15_rep1 = AUDPC.cal(reading.dates_15,wf_15_rep1)
AUDPC_wf_15_rep2 = AUDPC.cal(reading.dates_15,wf_15_rep2)


colnames(AUDPC_b_14_rep1) = c('ind','AUDPC_b_14_rep1')
colnames(AUDPC_b_14_rep2) = c('ind','AUDPC_b_14_rep2')
colnames(AUDPC_bf_14_rep1) = c('ind','AUDPC_bf_14_rep1')
colnames(AUDPC_bf_14_rep2) = c('ind','AUDPC_bf_14_rep2')
colnames(AUDPC_bf_14_rep3) = c('ind','AUDPC_bf_14_rep3')
colnames(AUDPC_bf_14_rep4) = c('ind','AUDPC_bf_14_rep4')
colnames(AUDPC_wf_14_rep1) = c('ind','AUDPC_wf_14_rep1')
colnames(AUDPC_wf_14_rep2) = c('ind','AUDPC_wf_14_rep2')
colnames(AUDPC_wf_14_rep3) = c('ind','AUDPC_wf_14_rep3')
colnames(AUDPC_wf_14_rep4) = c('ind','AUDPC_wf_14_rep4')
colnames(AUDPC_b_15_rep1) = c('ind','AUDPC_b_15_rep1')
colnames(AUDPC_b_15_rep2) = c('ind','AUDPC_b_15_rep2')
colnames(AUDPC_bf_15_rep1) = c('ind','AUDPC_bf_15_rep1')
colnames(AUDPC_bf_15_rep2) = c('ind','AUDPC_bf_15_rep2')
colnames(AUDPC_wf_15_rep1) = c('ind','AUDPC_wf_15_rep1')
colnames(AUDPC_wf_15_rep2) = c('ind','AUDPC_wf_15_rep2')

test = cbind(AUDPC_b_14_rep1,
             AUDPC_b_14_rep2,
             AUDPC_bf_14_rep1,
             AUDPC_bf_14_rep2,
             AUDPC_bf_14_rep3,
             AUDPC_bf_14_rep4,
             AUDPC_wf_14_rep1,
             AUDPC_wf_14_rep2,
             AUDPC_wf_14_rep3,
             AUDPC_wf_14_rep4,
             AUDPC_b_15_rep1,
             AUDPC_b_15_rep2,
             AUDPC_bf_15_rep1,
             AUDPC_bf_15_rep2,
             AUDPC_wf_15_rep1,
             AUDPC_wf_15_rep2)
 
 





