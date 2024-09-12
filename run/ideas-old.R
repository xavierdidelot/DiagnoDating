rm(list=ls())

set.seed(0)
library(ValidateDating,quietly=T)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10)

#perfect residuals
r0=resDating(dt,phy,resample=F,rate=10)
plotResid(r0);title(testResid(r0)$p.value)

#perform BactDating run
r=runDating(phy,dates,showProgress=T,nbIts=2e4)

#using bactdate run to resample
rtree=r$tree
rtree$root.time=NULL
k=sum(phy$edge.length)/sum(rtree$edge.length)
rtree$edge.length=rtree$edge.length*k
roottotip(rtree,dates)
r2=runDating(rtree,dates,minbralen=1e-10,model='strictgamma',showProgress=T,initMu=k,updateMu=F,updateRoot=F,useCoalPrior=T)#,initAlpha=mean(r$record[501:1000,'alpha']),updateAlpha=F
tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')
r3=resDating(takeSample(r2)$tree,phy,model='poisson',rate=r$rate,resample=F)
plotResid(r3);title(testResid(r3)$p.value)
