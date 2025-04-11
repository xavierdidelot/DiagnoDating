rm(list=ls())
library(ValidateDating,quietly=T)

set.seed(1)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

keepRoot=F;phy=unroot(phy)
rate=NA

r1=runDating(phy,dates,keepRoot=keepRoot,rate=rate)
v1=validate(r1);p1=median(v1$ps)
v1r=validate(r1,resampling = 2);p1r=median(v1r$ps)

r2=runDating(phy,dates,algo='treedater',keepRoot=keepRoot,rate=rate)
v2=validate(r2);p2=median(v2$ps)

r3=runDating(phy,dates,algo='node.dating',keepRoot=keepRoot,rate=rate)
v3=validate(r3);p3=median(v3$ps)

r4=runDating(phy,dates,algo='TreeTime',keepRoot=keepRoot,rate=rate)
set.seed(0);v4=validate(r4);p4=median(v4$ps)

r5=runDating(phy,dates,algo='LSD',keepRoot=keepRoot,rate=rate)
v5=validate(r5);p5=median(v5$ps)
