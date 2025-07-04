rm(list=ls())

library(DiagnoDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)
algos=c('BactDating','treedater','node.dating','TreeTime','LSD')

reps=100
allres <- foreach (rep = 1:(4*reps),.packages = c('ape','BactDating','DiagnoDating')) %dopar% {
  seed=floor(seq(1,reps+0.75,0.25))[rep]
  set.seed(seed)
  dates=runif(200,2000,2020)
  dt=simcoaltree(dates,10)
  phy=simobsphy(dt,mu=10,model='poisson')
  r0=resDating(dt,phy)
  p0=testResid(r0)$p.value
  if (rep%%4 == 1 || rep%%4 == 2) keepRoot=T else {keepRoot=F;phy=unroot(phy)}
  if (rep%%4 == 1 || rep%%4 == 3) rate=NA else rate=10

  l=list(seed=seed,keepRoot=keepRoot,rate=rate,truth=p0)

  for (j in 1:length(algos)) {
    r1=runDating(phy,dates,algo=algos[j],keepRoot=keepRoot,rate=rate)
    p1p=ppcheck(r1)
    if (j>1) r1=resample(r1)
    v1=postdistpvals(r1);p1=median(v1$ps)
    l=c(l,list(p1p,p1))
  }
  return(l)
}

parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
nam=names(allres[[1]])
allres=t(matrix(unlist(allres),nrow=14))
colnames(allres)=nam

w=which(allres[,'keepRoot']==F&is.na(allres[,'rate']))
hist(allres[w,6],main='',xlab='P-values',breaks=seq(0,1,0.05))

tab=matrix(NA,10,4)
for (i in 1:10) for (j in 1:4) {
  if (j==1 || j==2) keepRoot=F else keepRoot=T
  if (j==1 || j==3) isna=T else isna=F
  w=which(allres[,'keepRoot']==keepRoot&is.na(allres[,'rate'])==isna)
  tab[i,j]=length(which(allres[w,i+4]<0.05))
}
a2=rep(algos,each=2);a2[2*(1:length(algos))]=''
tab=cbind(a2,rep(c('PPcheck','Residuals'),5),tab)
tab=rbind(c('Method','Test','Given nothing','Given root','Given rate','Given root and rate'),tab)

library(xtable)
xt=xtable(tab,caption='Number of false positives found amongst a set of 100 replicates. Five different methods were used (BactDating, treedater, node.dating, TreeTime and LSD)
and two different tests (posterior predictive check and residual analysis). Each method was applied in four different conditions: given the correct root, given the
correct rate, given both or given neither.',label='tab:falsePos')
align(xt)<-'|r|r|r|r|r|r|r|'
print(xt,file='falsePos.tex',include.rownames=FALSE,include.colnames=FALSE,hline.after = c(0,1,nrow(xt)))
