rm(list=ls())

set.seed(0)
library(DiagnoDating,quietly=T)
dates=runif(50,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

#perfect residuals
r0=resDating(dt,phy,rate=10)
plotResid(r0);title(testResid(r0)$p.value)#great

#perform dating
r=runDating(phy,dates,algo='node.dating')
plotResid(r);title(testResid(r)$p.value)#underdispersed

#compute posterior distribution of p-values
rr=resample(r)
p=postdistpvals(rr,showPlot = T)
plotResid(p$last);title(testResid(p$last)$p.value)
