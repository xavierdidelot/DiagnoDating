library(DiagnoDating,quietly=T)
library(phangorn)
rm(list=ls())
source('fig_label.R')
source('roottotip.R')
set.seed(11)
dates=runif(100,2010,2020)
dt=simcoaltree(dates,alpha=1)
tree=simobsphy(dt,mu=10,model='poisson')
w=sample(1:100,size=5,replace=F)
w2=which(tree$edge[,2]%in%w)
tree$edge.length[w2]=tree$edge.length[w2]+20#runif(5,20,20)
tree2=drop.tip(tree,w)
#r=runDating(tree,dates,algo='BactDating',model='poisson',showProgress=T,keepRoot = T)

library(treedater)
tre=tree
tre$edge.length=tre$edge.length/1000
sts=dates
names(sts)=tre$tip.label
td=dater(tre,sts)
o=outlierLineages(td,alpha=0.01)
#goodnessOfFitPlot(td)
r=resDating(td,tree)

tre=tree2
tre$edge.length=tre$edge.length/1000
sts=dates
names(sts)=tree$tip.label
td2=dater(tre,sts)
o2=outlierLineages(td2,alpha=0.01)
#goodnessOfFitPlot(td)
r2=resDating(td2,tree2)


pdf('outlier.pdf',8,8)
par(mfrow=c(2,2),mar=c(4,4,2,2))
z=roottotip(tree,dates,rate=10)
fig_label('B',cex=2)
plotLikBranches(r,sub=1,minProb=1e-3)
fig_label('C',cex=2)
plotLikBranches(r2,sub=1,minProb=1e-3)
fig_label('D',cex=2)
dev.off()
system('open outlier.pdf')
