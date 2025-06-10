rm(list=ls())
#set.seed(0)
mu=2#known
alpha=10#known
dates=runif(300,2000,2020)
ttree=BactDating::simcoaltree(dates,alpha=alpha)
l=ttree$edge.length
n=length(l)
s=rpois(n,l*mu)
par(mfrow=c(2,3))

#u=runif(n,ppois(s-1,l*mu),ppois(s,l*mu))
#r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
#hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Exact residuals p=%.3g',r$p.value))

phy=ttree
phy$edge.length=s
#infer=runDating(phy,dates,initMu=mu,updateMu=F,initAlpha=alpha,updateAlpha=F,model='poisson',updateRoot=F,showProgress=T);infer=DiagnoDating:::reorderEdges(infer$tree,phy)
d<-ape::estimate.dates(phy,dates,mu=mu,lik.tol=1e-5);infer=phy;for (i in 1:nrow(infer$edge)) infer$edge.length[i]=d[infer$edge[i,2]]-d[infer$edge[i,1]]
l2=infer$edge.length

u=runif(n,ppois(s-1,l2*mu),ppois(s,l2*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('No resampling p=%.3g',r$p.value))

m=mean(l2);v=var(l2)
#k=m*m/v;theta=v/m#simpler guessing prior
k=m^2*mu/(v*mu-m);theta=v/m-1/mu#guessing prior using the Law of total expectation and the Law of total variance
l3=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l3*mu),ppois(s,l3*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('IID residuals p=%.3g',r$p.value))

if (T) {
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
}

s2=rpois(n,l2*mu)
phy2=ttree
phy2$edge.length=s2
d<-ape::estimate.dates(phy2,dates,mu=mu,lik.tol=1e-5);infer=phy2;for (i in 1:nrow(infer$edge)) infer$edge.length[i]=d[infer$edge[i,2]]-d[infer$edge[i,1]]
l5=infer$edge.length
u=runif(n,ppois(s-1,l5*mu),ppois(s,l5*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Resampling with parametric bootstrap p=%.3g',r$p.value))

l6=l5#(l5*4+l2)/5
u=runif(n,ppois(s-1,l6*mu),ppois(s,l6*mu))
u=u[1-u>1e-5]#remove extreme values
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Resampling with partial parametric bootstrap p=%.3g',r$p.value))

l7=(l5*4+l2)/5
u=runif(n,ppois(s-1,l7*mu),ppois(s,l7*mu))
r=DescTools::AndersonDarlingTest(u,null='punif',min=0,max=1)
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main=sprintf('Resampling with partial parametric bootstrap p=%.3g',r$p.value))
