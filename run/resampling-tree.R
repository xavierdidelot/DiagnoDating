rm(list=ls())
#set.seed(0)
mu=10#known
alpha=10#known
dates=runif(300,2000,2020)
ttree=BactDating::simcoaltree(dates,alpha=alpha)
l=ttree$edge.length
n=length(l)
s=rpois(n,l*mu)
par(mfrow=c(2,2))

u=runif(n,ppois(s-1,l*mu),ppois(s,l*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Exact residuals p=%.3g',r$p.value))

phy=ttree
phy$edge.length=s
#infer=runDating(phy,dates,initMu=mu,updateMu=F,initAlpha=alpha,updateAlpha=F,model='poisson',updateRoot=F,showProgress=T);infer=ValidateDating:::reorderEdges(infer$tree,phy)
d<-ape::estimate.dates(phy,dates,mu=mu);infer=phy;for (i in 1:nrow(infer$edge)) infer$edge.length[i]=d[infer$edge[i,2]]-d[infer$edge[i,1]]
l2=infer$edge.length

#mle=s/mu;v=0.01;l2=rgamma(n,shape=mle^2/v,scale=v/mle)#point estimate of l, centered on MLE
u=runif(n,ppois(s-1,l2*mu),ppois(s,l2*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('No resampling p=%.3g',r$p.value))

m=mean(l2);v=var(l2)
k=m*m/v;theta=v/m#simpler guessing prior
l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('IID residuals p=%.3g',r$p.value))

store=matrix(NA,100,n)
for (i in 1:nrow(store)) {
  t=BactDating::simcoaltree(dates,alpha=alpha)$edge.length
  s2=rpois(n,t*mu)
  ind=order(s2+runif(length(s2),0,0.01))
  t=t[ind]
  s2=s2[ind]
  store[i,]=t
}

m=v=k=theta=rep(NA,n)
ind=order(l2+runif(length(l2),0,1e-5))
ind[ind]=1:n
for (j in 1:n) {
  m[j]=mean(as.vector(store[,ind[j]]))
  v[j]=var (as.vector(store[,ind[j]]))
  k[j]=m[j]*m[j]/v[j]
  theta[j]=v[j]/m[j]
}

l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Coalescent resampling p=%.3g',r$p.value))

