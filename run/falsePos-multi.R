rm(list=ls())

library(ValidateDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

reps=100
allres <- foreach (rep = 1:(4*reps),.packages = c('ape','BactDating','ValidateDating')) %dopar% {
  seed=floor(seq(1,reps+0.75,0.25))[rep]
  set.seed(seed)
  dates=runif(200,2000,2020)
  dt=simcoaltree(dates,10)
  phy=simobsphy(dt,mu=10,model='poisson')
  r0=resDating(dt,phy)
  p0=testResid(r0)$p.value
  if (rep%%4 == 1 || rep%%4 == 2) keepRoot=T else {keepRoot=F;phy=unroot(phy)}
  if (rep%%4 == 1 || rep%%4 == 3) rate=NA else rate=10
  r1=runDating(phy,dates,keepRoot=keepRoot,rate=rate)
  v1=validate(r1,resampling = 0);p1=median(v1$ps)
  v1r=validate(r1,resampling = 2);p1r=median(v1r$ps)
  r2=runDating(phy,dates,algo='treedater',keepRoot=keepRoot,rate=rate)
  v2=validate(r2,resampling = 2);p2=median(v2$ps)
  r3=runDating(phy,dates,algo='node.dating',keepRoot=keepRoot,rate=rate)
  v3=validate(r3,resampling = 2);p3=median(v3$ps)
  r4=runDating(phy,dates,algo='TreeTime',keepRoot=keepRoot,rate=rate)
  v4=validate(r4,resampling = 2);p4=median(v4$ps)
  return(list(seed=seed,keepRoot=keepRoot,rate=rate,truth=p0,bactdate=p1,bactdate2=p1r,treedater=p2,nodedating=p3,treetime=p4))
}

parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
nam=names(allres[[1]])
allres=t(matrix(unlist(allres),nrow=9))
colnames(allres)=nam

w=which(allres[,'keepRoot']==F&!is.na(allres[,'rate']))
hist(allres[w,'truth'],main='',xlab='P-values',breaks=seq(0,1,0.1))
