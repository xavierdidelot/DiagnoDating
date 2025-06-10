rm(list=ls())

set.seed(0)
library(DiagnoDating,quietly=T)
dates=runif(50,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=1,model='poisson')

#perfect residuals
r0=resDating(dt,phy,rate=1)
plotResid(r0);title(testResid(r0)$p.value)#great

#perform dating
r=runDating(phy,dates,algo='node.dating')
plotResid(r);title(testResid(r)$p.value)#underdispersed

#compute posterior distribution of p-values
validate(r,resampling = 2,showLast=T)

#compute posterior distribution of p-values
validate(r,resampling = 3,showLast=T)
