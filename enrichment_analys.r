all_gene = read.csv('all_gene.csv', header = F)
rgene = read.csv('R_gene.csv', header = F)
qtl_gene = read.csv('QTL_gene.csv', header = F)
qtl_r_gene = read.csv("QTL_R_gene.csv", header = F)
N = length(all_gene[,1])    #number of all gene
n = length(qtl_gene[,1])    #number of genes in QTL  (genes sampled)
M_all = length(rgene[,1])       #number of resistance related gene (all classes)
k_all = length(qtl_r_gene[,1])  #number of resistance related gene in QTL (all classes)

gbyf = split(rgene,rgene$V2)  #split differen Rgene by function
gbyf_qtl = split(qtl_r_gene,qtl_r_gene$V2)
res = data.frame(matrix(NA,nrow = 9,ncol = 6))
names(res)=c("Gene class", "Pvalue", "actual/expected","expected Rgene","Rgene in QTL", "All Rgene")

index = 1

for (i in names(gbyf)) {
  #list = gbyf[[i]]
  #browser()
  M = length(gbyf[[i]][,2])  #total R gene
  expd_count = n*(M/N)     #expected Rgene of one class in QTL region, if sample randomed
  k = length(gbyf_qtl[[i]][,2]) #number of R gene in QTL
  OddsRatio = k/expd_count # actual Rgene / expected Rgene in QTL regions
  p = phyper(k, M, N-M, n,lower.tail = F)
  print(paste(i,p,OddsRatio,expd_count,k,M,sep="    "))
  res[index,]=c(i,p,OddsRatio,expd_count,k,M)
  index = index + 1
}

write.csv(res,'hypogeo.test.result.csv')
