rm(list=ls())

set.seed(4)
library(ValidateDating,quietly=T)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

keepRoot=T
rate=10

#r2=runDating(phy,dates,algo='treedater',keepRoot=keepRoot,rate=rate)
#v2=validate(r2,nbIts=1000);p2=median(v2$ps)

#r3=runDating(phy,dates,algo='node.dating',keepRoot=keepRoot,rate=rate)
#v3=validate(r3);p3=median(v3$ps)

r4=runDating(phy,dates,algo='TreeTime',keepRoot=keepRoot,rate=rate)
set.seed(0);v4=validate(r4);p4=median(v4$ps)
