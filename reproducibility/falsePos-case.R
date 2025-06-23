rm(list=ls())
library(DiagnoDating,quietly=T)

set.seed(1)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

keepRoot=F;phy=unroot(phy)
rate=NA

r1=runDating(phy,dates,keepRoot=keepRoot,rate=rate)
p1p=ppcheck(r1)
v1=postdistpvals(r1);p1=median(v1$ps)

r2=runDating(phy,dates,algo='treedater',keepRoot=keepRoot,rate=rate)
p2p=ppcheck(r2)
r2=resample(r2);v2=postdistpvals(r2);p2=median(v2$ps)

r3=runDating(phy,dates,algo='node.dating',keepRoot=keepRoot,rate=rate)
p3p=ppcheck(r3)
r3=resample(r3);v3=postdistpvals(r3);p3=median(v3$ps)

r4=runDating(phy,dates,algo='TreeTime',keepRoot=keepRoot,rate=rate)
p4p=ppcheck(r4)
r4=resample(r4);v4=postdistpvals(r4);p4=median(v4$ps)

r5=runDating(phy,dates,algo='LSD',keepRoot=keepRoot,rate=rate)
p5p=ppcheck(r5)
r5=resample(r5);v5=postdistpvals(r5);p5=median(v5$ps)
