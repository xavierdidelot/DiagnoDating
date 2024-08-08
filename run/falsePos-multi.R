rm(list=ls())

library(ValidateDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

allres <- foreach (rep = 1:100,.packages = c('ape','BactDating','ValidateDating')) %dopar% {
  set.seed(rep)
  dates=runif(200,2000,2020)
  dt=simcoaltree(dates,10)
  phy=simobsphy(dt,mu=10)
  r0=resDating(dt,phy,resample=F)
  p0=testResid(r0)$p.value
  r1=runDating(phy,dates)
  p1=testResid(r1)$p.value
  r2=runDating(phy,dates,algo='treedater')
  p2=testResid(r2)$p.value
  r3=runDating(phy,dates,algo='node.dating')
  p3=testResid(r3)$p.value
  r4=runDating(phy,dates,algo='TreeTime')
  p4=testResid(r4)$p.value
  return(list(seed=rep,p0,p1,p2,p3,p4))
}
parallel::stopCluster(cl)
save.image('falsePos-multi.RData')
