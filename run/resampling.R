rm(list=ls())
mu=2#known
n=1e6
shape=2
scale=2
l=rgamma(n,shape=shape,scale=scale)
s=rpois(length(l),l*mu)
par(mfrow=c(2,2))

u=runif(n,ppois(s-1,l*mu),ppois(s,l*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Exact residuals')

mle=s/mu;v=0.1;l2=rgamma(n,shape=mle^2/v,scale=v/mle)#point estimate of l, centered on MLE
u=runif(n,ppois(s-1,l2*mu),ppois(s,l2*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a point estimate without resampling')

l3=rgamma(n,shape=shape+s,scale=scale/(1+scale*mu))#sample from posterior (or could use multiple and use posterior distribution of p-values)
u=runif(n,ppois(s-1,l3*mu),ppois(s,l3*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a posterior sample')

#k=shape;theta=scale#prior known
m=mean(mle);v=var(mle);k=m^2/v;theta=v/m#guessing prior
l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a point estimate with resampling')
