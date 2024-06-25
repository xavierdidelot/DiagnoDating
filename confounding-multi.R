library(BactDating)
library(ape)
library(mlesky)
source('simStructure.R')

allres=matrix(NA,10,3)
for (rep in 1:10) {
  set.seed(rep)
  s=list(runif(30,2010,2010.01),runif(30,2010.1,2010.11),runif(30,2010.2,2010.21))
  dt=simStructure(popStarts=c(2000,2000,2000),globalNeg = 30,samplingDates=s)
  dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

  tree=simobsphy(dt,mu=10,model='poisson')
  res=suppressWarnings(roottotip(tree,dates,showFig = F))
  allres[rep,]=c(res$rate,res$ori-dt$root.time,res$pvalue)
}
w=which(allres[,3]<0.01)
par(mfrow=c(2,1))
hist(allres[w,1],xlab='Clock rate',main='')
hist(allres[w,2],xlab='Root date',main='')

