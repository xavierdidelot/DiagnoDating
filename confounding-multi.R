library(BactDating)
library(ape)
source('simStructure.R')
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating')) %dopar% {
  set.seed(rep)
  s=list(runif(30,2010,2010.01),runif(30,2010.1,2010.11),runif(30,2010.2,2010.21))
  dt=simStructure(popStarts=c(2000,2000,2000),globalNeg = 10,samplingDates=s)
  dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

  tree=simobsphy(dt,mu=10,model='poisson')
  res=suppressWarnings(roottotip(tree,dates,showFig = F))
  r=bactdate(tree,dates,model='poisson',updateRoot = F)
  return(list(seed=rep,rate=res$rate,ori=res$ori-dt$root.time,pval=res$pvalue,resbactdate=r))
}
parallel::stopCluster(cl)
save.image('confounding-multi.RData')

v1=c()
v2=c()
for (i in 1:length(allres)) {
  if (allres[[i]]$pval<0.01) {
    v1=c(v1,allres[[i]]$rate)
    v2=c(v2,allres[[i]]$ori)
  }
}
par(mfrow=c(2,1))
hist(v1,xlab='Clock rate (truth is 10)',main='',breaks=20) 
hist(v2,xlab='Root date (truth is 0)',main='',breaks=20)

