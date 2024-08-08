rm(list=ls())

set.seed(10)
library(ValidateDating,quietly=T)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10)
r0=resDating(dt,phy)
print(testResid(r0))
r=runDating(phy,dates,showProgress=T,nbIts=2e4)
plotResid(r)
print(testResid(r))

r2=runDating(phy,dates,algo='treedater')
plotResid(r2)
print(testResid(r2))

r3=runDating(phy,dates,algo='node.dating')
plotResid(r3)
print(testResid(r3))

r4=runDating(phy,dates,algo='TreeTime')
plotResid(r4)
print(testResid(r4))

