rm(list=ls())

set.seed(5)
library(ValidateDating,quietly=T)
dates=runif(200,2000,2020)
dt=simcoaltree(dates,10)
phy=simobsphy(dt,mu=10,model='poisson')

r1=runDating(phy,dates,algo='node.dating')
validate(r1,resampling=2,showLast = T)

r2=runDating(phy,dates,algo='node.dating',rate=10)
validate(r2,resampling=2,showLast = T)


