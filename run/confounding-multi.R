library(BactDating)
library(ape)
source('simStructure.R')
library(foreach)
source('wrappers.R')

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating')) %dopar% {
  set.seed(rep)
  s=list(runif(30,2010,2010.01),runif(30,2010.1,2010.11),runif(30,2010.2,2010.21))
  dt=simStructure(popStarts=c(2000,2000,2000),globalNeg = 10,samplingDates=s)
  dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

  tree=simobsphy(dt,mu=10,model='poisson')
  res=suppressWarnings(roottotip(tree,dates,showFig = F))
  rbactdate=bactdate(tree,dates,model='poisson',updateRoot = F)
  rlsd=runLSD(tree,dates,rep)
  return(list(seed=rep,dt=dt,rate=res$rate,ori=res$ori-dt$root.time,pval=res$pvalue,resbactdate=rbactdate,lsdrate=rlsd$rate,lsdtmrca=rlsd$tmrca))
}
parallel::stopCluster(cl)
save.image('confounding-multi.RData')

rm(list=ls())
load('confounding-multi.RData')
v1=v2=v3=v4=v5=v6=c()
for (i in 1:length(allres)) {
  if (allres[[i]]$pval<0.01) {
    v1=c(v1,allres[[i]]$rate)
    v2=c(v2,allres[[i]]$ori)
  }
  rec=allres[[i]]$resbactdate$record
  rec=rec[(nrow(rec)/2):nrow(rec),]
  v3=c(v3,mean(rec[,'mu']))
  v4=c(v4,allres[[i]]$resbactdate$rootdate-allres[[i]]$dt$root.time)
  v5=c(v5,allres[[i]]$lsdrate)
  v6=c(v6,allres[[i]]$lsdtmrca-allres[[i]]$dt$root.time)
}
par(mfrow=c(3,2))
hist(v1,xlab='Clock rate (truth is 10)',main='',breaks=20) 
hist(v2,xlab='Root date (truth is 0)',main='',breaks=20)
hist(v3,xlab='Clock rate (truth is 10)',main='',breaks=20) 
hist(v4,xlab='Root date (truth is 0)',main='',breaks=20)
hist(v5,xlab='Clock rate (truth is 10)',main='',breaks=20) 
hist(v6,xlab='Root date (truth is 0)',main='',breaks=20)
