rm(list=ls())

library(DiagnoDating,quietly=T)
algos=c('BactDating','treedater','node.dating','TreeTime','LSD')

set.seed(1)#
popStarts=c(2009,2010,2011)-50
v1=60;v2=60.1
s=list(runif(50,popStarts[1]+v1,popStarts[1]+v1),runif(50,popStarts[2]+v1,popStarts[2]+v2),runif(50,popStarts[3]+v1,popStarts[3]+v2))
dt=simStructure(popStarts = popStarts,samplingDates=s,globalNeg = 10)
dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

tree=simobsphy(dt,mu=10,model='poisson')
#tree=unroot(tree)
rd=runDating(tree,dates,showProgress=T,nbIts=2e4,initMu=15)
rd$p1=ppcheck(rd)
rd$p2=median(postdistpvals(rd)$ps)

bd=rd
class(bd)<-'resBactDating'
plot(bd,'trace')

roottotip(initRoot(tree,dates),dates)

source('fig_label.R')
source('roottotip.R')
pdf('confounding-case.pdf',8,8)
par(mfrow=c(2,2),mar=c(4,4,2,2))
z=roottotip(initRoot(tree,dates),dates)
fig_label('B',cex=2)
#plotLikBranches(rd,sub=1,minProb=1e-3)
plotResid(postdistpvals(rd)$last,4,main='',xlim=c(-3,3),ylim=c(-4,4))
fig_label('C',cex=2)
ps=postdistpvals(rd)$ps
h=hist(ps,breaks=seq(0,1,length.out=51),xlab='',ylab='',main='')
y=length(which(ps<0.05))
par(xpd=NA)
#text(0.025,y,sprintf('%.1f%%',y*100/length(ps)),pos=3)
med=median(ps)
lines(c(med,med),c(0,max(h$counts)),lty=2)
#text(med,max(h$counts),format(med,digits=3),pos=4)
par(xpd=F)
fig_label('D',cex=2)
dev.off()
system('open confounding-case.pdf')
