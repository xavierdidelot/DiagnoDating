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
  p0=length(which(validate(r0)<0.05))
  r1=runDating(phy,dates)
  p1=length(which(validate(r1)<0.05))
  r2=runDating(phy,dates,algo='treedater')
  p2=length(which(validate(r2)<0.05))
  r3=runDating(phy,dates,algo='node.dating')
  p3=length(which(validate(r3)<0.05))
  r4=runDating(phy,dates,algo='TreeTime')
  p4=length(which(validate(r4)<0.05))
  return(list(seed=rep,p0,p1,p2,p3,p4))
}
parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
allres=t(matrix(unlist(allres),nrow=6))
