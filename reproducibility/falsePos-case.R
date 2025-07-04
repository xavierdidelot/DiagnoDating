rm(list=ls())
library(DiagnoDating,quietly=T)
algos=c('BactDating','treedater','node.dating','TreeTime','LSD')

set.seed(1)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

keepRoot=F;phy=unroot(phy)
rate=NA

lr=list()
lpp=list()
lp=list()
for (j in 1:length(algos)) {
  r1=runDating(phy,dates,algo=algos[j],keepRoot=keepRoot,rate=rate)
  p1p=ppcheck(r1)
  if (j>1) r1=resample(r1)
  v1=postdistpvals(r1);p1=median(v1$ps)
  lr=append(lr,r1)
  lpp=append(lpp,p1p)
  lp=append(lp,p1)
}
