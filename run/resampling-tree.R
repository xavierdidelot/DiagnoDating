rm(list=ls())
mu=10#known
alpha=10
dates=runif(500,2000,2020)
ttree=BactDating::simcoaltree(dates,alpha=alpha)
l=ttree$edge.length
n=length(l)
s=rpois(n,l*mu)
par(mfrow=c(2,2))

u=runif(n,ppois(s-1,l*mu),ppois(s,l*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Exact residuals')
DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)

phy=ttree
phy$edge.length=s
#infer=runDating(phy,dates,initMu=mu,updateMu=F,initAlpha=alpha,updateAlpha=F,model='poisson',updateRoot=F,showProgress=T);infer=ValidateDating:::reorderEdges(infer$tree,phy)
d<-ape::estimate.dates(phy,dates,mu=mu);infer=phy;for (i in 1:nrow(infer$edge)) infer$edge.length[i]=d[infer$edge[i,2]]-d[infer$edge[i,1]]
l2=infer$edge.length

#mle=s/mu;v=0.01;l2=rgamma(n,shape=mle^2/v,scale=v/mle)#point estimate of l, centered on MLE
u=runif(n,ppois(s-1,l2*mu),ppois(s,l2*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='No resampling')
DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)

m=mean(l2);v=var(l2)
k=m*m/v;theta=v/m#simpler guessing prior
l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='IID resampling')
DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)

store=matrix(NA,100,n)
for (i in 1:nrow(store)) {
  t=BactDating::simcoaltree(dates,alpha=alpha)$edge.length
  s2=rpois(n,t*mu)
  ind=order(s2)
  t=t[ind]
  store[i,]=t
}

m=v=k=theta=rep(NA,n)
ind=order(l2)
ind[ind]=1:n
for (j in 1:n) {
  m[j]=mean(as.vector(store[,ind[j]]))
  v[j]=var (as.vector(store[,ind[j]]))
  k[j]=m[j]*m[j]/v[j]
  theta[j]=v[j]/m[j]
  #k[j]=m[j]^2*mu/(v[j]*mu-m[j])
  #theta[j]=v[j]/m[j]-1/mu
}

l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Real resampling')
DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
