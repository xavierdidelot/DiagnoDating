rm(list=ls())

library(DiagnoDating,quietly=T)
algos=c('BactDating','treedater','node.dating','TreeTime','LSD')

set.seed(14)#also 14 100
popStarts=c(2009,2010,2011)-10
v1=20;v2=20.1
s=list(runif(100,popStarts[1]+v1,popStarts[1]+v1),runif(100,popStarts[2]+v1,popStarts[2]+v2),runif(100,popStarts[3]+v1,popStarts[3]+v2))
dt=simStructure(popStarts = popStarts,samplingDates=s,globalNeg = 10)
dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

tree=simobsphy(dt,mu=10,model='poisson')
tree=unroot(tree)
l=list(seed=rep,dt=dt)
for (i in 1:5) {
    rd=runDating(tree,dates,algo=algos[i])
    if (i==1) {
      rd$p1=ppcheck(rd)
      rd$p2=median(postdistpvals(rd)$ps)
    }
  l=c(l,list(rd))
}

bd=l[[3]]
class(bd)<-'resBactDating'
plot(bd,'trace')

roottotip(initRoot(tree,dates),dates)
