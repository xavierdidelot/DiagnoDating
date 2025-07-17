rm(list=ls())

library(DiagnoDating,quietly=T)
library(doParallel)

cl <- makeCluster(parallel::detectCores() - 1, type = "PSOCK", outfile="")
registerDoParallel(cl)
algos=c('BactDating','treedater','node.dating','TreeTime','LSD')

reps=100
allres <- foreach (rep = 1:reps,.packages = c('ape','BactDating','DiagnoDating')) %dopar% {
  cat('Starting task',rep,'\n')
  capture.output({
    set.seed(rep)
    popStarts=c(2009,2010,2011)-10
    v1=20;v2=20.1
    s=list(runif(100,popStarts[1]+v1,popStarts[1]+v1),runif(100,popStarts[2]+v1,popStarts[2]+v2),runif(100,popStarts[3]+v1,popStarts[3]+v2))
    dt=simStructure(popStarts = popStarts,samplingDates=s,globalNeg = 10)
    dates=unname(dist.nodes(dt)[Ntip(dt)+1,1:Ntip(dt)])+dt$root.time

    tree=simobsphy(dt,mu=10,model='poisson')
    tree=unroot(tree)
    l=list(seed=rep,dt=dt)
    for (i in 1:5) {
      if (i==1){
      rd=runDating(tree,dates,algo=algos[i])
      if (i==1) {
        rd$p1=ppcheck(rd)
        rd$p2=median(postdistpvals(rd)$ps)
      }}
      l=c(l,list(rd))
    }
  })
  cat('Finished task',rep,'\n')
  return(l)
}
parallel::stopCluster(cl)
save.image('confounding-multi.RData')

rm(list=ls())
load('confounding-multi.RData')
v=matrix(NA,length(allres),10)
for (i in 1:length(allres)) {
  for (j in 1:5) {
    res=allres[[i]][[j+2]]
    v[i,1+(j-1)*2  ]=res$rate
    v[i,1+(j-1)*2+1]=res$rootdate-allres[[i]]$dt$root.time
  }
}
colnames(v)=sprintf(rep(c('Rate %s','TMRCA %s'),5),rep(algos,each=2))
pdf('confounding-multi.pdf',10,10)
par(mfrow=c(5,2))
for (i in 1:10) {
  hist(v[,i],breaks=20,xlab = '',ylab='',main=colnames(v)[i])
}
dev.off()
system('open confounding-multi.pdf')

allpvals=matrix(NA,length(allres),2)
for (i in 1:length(allres)) {
  allpvals[i,1]=allres[[i]][[3]]$p1
  allpvals[i,2]=allres[[i]][[3]]$p2
}
