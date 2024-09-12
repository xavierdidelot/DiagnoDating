rm(list=ls())

set.seed(10)
library(ValidateDating,quietly=T)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10)
r0=resDating(dt,phy)#perfect residuals
plotResid(r0);title(testResid(r0)$p.value)

r=runDating(phy,dates,showProgress=T,nbIts=2e4)
plotResid(r);title(testResid(r)$p.value)
validate(r)

r2=runDating(phy,dates,algo='treedater')
plotResid(r2);title(testResid(r2)$p.value)

r3=runDating(phy,dates,algo='node.dating')
plotResid(r3);title(testResid(r3)$p.value)


r4=runDating(phy,dates,algo='TreeTime')
plotResid(r4);title(testResid(r4)$p.value)


