rm(list=ls())

library(ValidateDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

reps=100
allres <- foreach (rep = 1:(2*reps),.packages = c('ape','BactDating','ValidateDating')) %dopar% {
  set.seed(c(1:reps,1:reps)[rep])
  dates=runif(200,2000,2020)
  dt=simcoaltree(dates,10)
  phy=simobsphy(dt,mu=10,model='poisson')
  r0=resDating(dt,phy)
  p0=testResid(r0)$p.value
  if (rep<=reps) keepRoot=T else {keepRoot=F;phy=unroot(phy)}
  r1=runDating(phy,dates,keepRoot=keepRoot)
  v1=validate(r1,resampling = 0);p1=median(v1$ps)
  v1r=validate(r1,resampling = 2);p1r=median(v1r$ps)
  r2=runDating(phy,dates,algo='treedater',keepRoot=keepRoot)
  v2=validate(r2,resampling = 2);p2=median(v2$ps)
  r3=runDating(phy,dates,algo='node.dating',keepRoot=keepRoot)
  v3=validate(r3,resampling = 2);p3=median(v3$ps)
  r4=runDating(phy,dates,algo='TreeTime',keepRoot=keepRoot)
  v4=validate(r4,resampling = 2);p4=median(v4$ps)
  r1=runDating(phy,dates,rate=10,keepRoot=keepRoot)
  v1=validate(r1,resampling = 0);p1f=median(v1$ps)
  v1r=validate(r1,resampling = 2);p1fr=median(v1r$ps)
  r2=runDating(phy,dates,algo='treedater',rate=10,keepRoot=keepRoot)
  v2=validate(r2,resampling = 2);p2f=median(v2$ps)
  r3=runDating(phy,dates,algo='node.dating',rate=10,keepRoot=keepRoot)
  v3=validate(r3,resampling = 2);p3f=median(v3$ps)
  r4=runDating(phy,dates,algo='TreeTime',rate=10,keepRoot=keepRoot)
  v4=validate(r4,resampling = 2);p4f=median(v4$ps)
  return(list(seed=rep,p0,p1,p1f,p2,p3,p4,p1f,p1fr,p2f,p3f,p4f))
  #                1.  2. 3. 4.  5. 6. 7. 8.  9.   10. 11. 12.
}

parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
allres=t(matrix(unlist(allres),nrow=12))
