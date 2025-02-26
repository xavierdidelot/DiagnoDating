rm(list=ls())

set.seed(0)
library(ValidateDating,quietly=T)
dates=runif(50,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

#perfect residuals
r0=resDating(dt,phy,rate=10)
plotResid(r0);title(testResid(r0)$p.value)#great

#perform dating
r=runDating(phy,dates,algo='node.dating')
plotResid(r);title(testResid(r)$p.value)#underdispersed

#using bactdate run to resample
rtree=r$tree
rtree$root.time=NULL
k=sum(phy$edge.length)/sum(rtree$edge.length)
rtree$edge.length=rtree$edge.length*k
roottotip(rtree,dates)#perfect
r2=runDating(rtree,dates,minbralen=1e-10,model='strictgamma',showProgress=T,initMu=k,updateMu=F,updateRoot=F)#,initAlpha=mean(r$record[501:1000,'alpha']),updateAlpha=F)
tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')#just to check traces

#check residuals for a single sample
nr=nrow(r2$record)
samtree=takeSample(r2,w=sample((nr/2):nr,1))$tree
r3=resDating(samtree,phy,algo=r$algo,model='poisson',rate=r$rate)
plotResid(r3);title(testResid(r3)$p.value)

#compute posterior distribution of p-values
validate(r,resampling = 2)
