rm(list=ls())

library(DiagnoDating,quietly=T)
library(foreach)

cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
doParallel::registerDoParallel(cl)

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

    return(list(seed=seed,keepRoot=keepRoot,rate=rate,truth=p0,bactdate0=p1p,bactdate=p1,treedater0=p2p,treedater=p2,nodedating0=p3p,nodedating=p3,treetime0=p4p,treetime=p4,lsd0=p5p,lsd=p5))
}

parallel::stopCluster(cl)
save.image('falsePos-multi.RData')

rm(list=ls())
load('falsePos-multi.RData')
nam=names(allres[[1]])
allres=t(matrix(unlist(allres),nrow=14))
colnames(allres)=nam

w=which(allres[,'keepRoot']==F&!is.na(allres[,'rate']))
hist(allres[w,'nodedating'],main='',xlab='P-values',breaks=seq(0,1,0.05))

tab=matrix(NA,10,4)
rownames(tab)=colnames(allres)[5:14]
colnames(tab)=c('Given nothing','Given root','Given rate','Given root and rate')
for (i in 1:10) for (j in 1:4) {
  if (j==1 || j==2) keepRoot=F else keepRoot=T
  if (j==1 || j==3) isna=T else isna=F
  w=which(allres[,'keepRoot']==keepRoot&is.na(allres[,'rate'])==isna)
  tab[i,j]=length(which(allres[w,i+4]<0.05))
}
