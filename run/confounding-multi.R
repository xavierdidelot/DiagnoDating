library(BactDating)
library(ape)
library(foreach)
library(confounding)

rm(list=ls())

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating','confounding')) %dopar% {
  set.seed(rep)
  s=list(runif(30,2010,2010.1),runif(30,2010.1,2010.2),runif(30,2010.2,2010.3))
  dt=simStructure(popStarts=c(2009,2009,2009),samplingDates=s)
  dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

  tree=simobsphy(dt,mu=10,model='poisson')
  r1=runDating(tree,dates,algo=1)
  r2=runDating(tree,dates,algo=2)
  r3=runDating(tree,dates,algo=3)
  r4=runDating(tree,dates,algo=4)
  return(list(seed=rep,dt=dt,r1,r2,r3,r4))
}
parallel::stopCluster(cl)
save.image('confounding-multi.RData')

rm(list=ls())
load('confounding-multi.RData')
v=matrix(NA,length(allres),8)
for (i in 1:length(allres)) {
    for (j in 1:4) {
      res=allres[[i]][[j+2]]
      v[i,1+(j-1)*2  ]=res$rate
      v[i,1+(j-1)*2+1]=res$rootdate-allres[[i]]$dt$root.time
    }
}
par(mfrow=c(4,2))
for (i in 1:8) {
  hist(v[,i],breaks=20)
}

