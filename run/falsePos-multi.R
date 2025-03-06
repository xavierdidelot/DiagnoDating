rm(list=ls())

library(ValidateDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating','ValidateDating')) %dopar% {
  set.seed(rep)
  dates=runif(200,2000,2020)
  dt=simcoaltree(dates,10)
  phy=simobsphy(dt,mu=10,model='poisson')
  r0=resDating(dt,phy)
  p0=testResid(r0)$p.value
  r1=runDating(phy,dates)
  v1=validate(r1,resampling = 0);p1=length(which(v1<0.05))/length(v1)
  r2=runDating(phy,dates,algo='treedater')
  v2=validate(r2,resampling = 2);p2=length(which(v2<0.05))/length(v2)
  r3=runDating(phy,dates,algo='node.dating')
  v3=validate(r3,resampling = 2);p3=length(which(v3<0.05))/length(v3)
  r4=runDating(phy,dates,algo='TreeTime')
  v4=validate(r4,resampling = 2);p4=length(which(v4<0.05))/length(v4)
  r1=runDating(phy,dates,rate=10)
  v1=validate(r1,resampling = 0);p1f=length(which(v1<0.05))/length(v1)
  r2=runDating(phy,dates,algo='treedater',rate=10)
  v2=validate(r2,resampling = 2);p2f=length(which(v2<0.05))/length(v2)
  r3=runDating(phy,dates,algo='node.dating',rate=10)
  v3=validate(r3,resampling = 2);p3f=length(which(v3<0.05))/length(v3)
  r4=runDating(phy,dates,algo='TreeTime',rate=10)
  v4=validate(r4,resampling = 2);p4f=length(which(v4<0.05))/length(v4)
  return(list(seed=rep,p0,p1,p2,p3,p4,p1f,p2f,p3f,p4f))
}
parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
allres=t(matrix(unlist(allres),nrow=6))
