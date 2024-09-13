rm(list=ls())

set.seed(0)
library(ValidateDating,quietly=T)
dates=runif(50,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

#perfect residuals
r0=resDating(dt,phy,rate=10)
plotResid(r0);title(testResid(r0)$p.value)

#perform dating
r=runDating(phy,dates,algo='node.dating')
r$inputtree=phy
plotResid(r);title(testResid(r)$p.value)

#using bactdate run to resample
rtree=r$tree
rtree$root.time=NULL
k=sum(phy$edge.length)/sum(rtree$edge.length)
rtree$edge.length=rtree$edge.length*k
roottotip(rtree,dates)#perfect
r2=runDating(rtree,dates,minbralen=1e-10,model='strictgamma',showProgress=T,initMu=k,updateMu=F,updateRoot=F)#,initAlpha=mean(r$record[501:1000,'alpha']),updateAlpha=F)
tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')#just to check traces
r3=resDating(r2$tree,phy,algo=r$algo,model='poisson',rate=r$rate)
r3$record=r2$record
r3$inputtree=r$inputtree
validate(r3,resampling = F)
