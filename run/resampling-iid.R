rm(list=ls())
mu=10#known
n=1e6
shape=2
scale=2
l=rgamma(n,shape=shape,scale=scale)
s=rpois(length(l),l*mu)
par(mfrow=c(2,3))

u=runif(n,ppois(s-1,l*mu),ppois(s,l*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Exact residuals')

l2=s/mu#point estimate of l, using MLE
if (F) {v=0.01;l2=rgamma(n,shape=l2^2/v,scale=v/l2)}#add a bit of randomness if desired
u=runif(n,ppois(s-1,l2*mu),ppois(s,l2*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a point estimate without resampling')

l3=rgamma(n,shape=shape+s,scale=scale/(1+scale*mu))#sample from posterior (or could use multiple and use posterior distribution of p-values)
u=runif(n,ppois(s-1,l3*mu),ppois(s,l3*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a posterior sample')

#k=shape;theta=scale#prior known
m=mean(l2);v=var(l2)
k=m^2*mu/(v*mu-m);theta=v/m-1/mu#guessing prior using the Law of total expectation and the Law of total variance
l4=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))#sample from posterior based on MLE
u=runif(n,ppois(s-1,l4*mu),ppois(s,l4*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a point estimate with resampling')

l5=suppressWarnings(rnorm(length(l2),l2,sqrt(l2^2/s)))
u=runif(n,ppois(s-1,l5*mu),ppois(s,l5*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a Laplace approximation')

s2=rpois(length(l),l2*mu)
l6=s2/mu
u=runif(n,ppois(s-1,l6*mu),ppois(s,l6*mu))
hist(u,breaks=seq(0,1,0.05),freq=F,xlab='',ylab='',main='Residuals based on a parametric bootstrap')
