library(BactDating)
library(ape)
library(foreach)
library(confounding)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating','confounding')) %dopar% {
  set.seed(rep)
  s=list(runif(30,2010,2010.1),runif(30,2010.1,2010.2),runif(30,2010.2,2010.3))
  dt=simStructure(popStarts=c(2009,2009,2009),samplingDates=s)
  dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

  tree=simobsphy(dt,mu=10,model='poisson')
  res=suppressWarnings(roottotip(tree,dates,showFig = F))
  rbactdate=bactdate(tree,dates,model='poisson',updateRoot = F)
  rlsd=runLSD(tree,dates)
  rtreedater=runTreeDater(tree,dates)
  return(list(seed=rep,dt=dt,rate=res$rate,ori=res$ori,pval=res$pvalue,rbactdate=rbactdate,rlsd=rlsd,rtreedater=rtreedater))
}
parallel::stopCluster(cl)
save.image('confounding-multi.RData')

rm(list=ls())
load('confounding-multi.RData')
v=matrix(NA,length(allres),8)
for (i in 1:length(allres)) {
    v[i,1]=allres[[i]]$rate
    v[i,2]=allres[[i]]$ori-allres[[i]]$dt$root.time
    for (j in 1:3) {
      rec=allres[[i]][[j+5]]$record
      rec=rec[(nrow(rec)/2):nrow(rec),]
      v[i,3+(j-1)*2. ]=mean(rec[,'mu'])
      v[i,3+(j-1)*2+1]=allres[[i]][[j+5]]$rootdate-allres[[i]]$dt$root.time
    }
}
par(mfrow=c(4,2))
for (i in 1:8) {
  hist(v[,i],breaks=20)
}
