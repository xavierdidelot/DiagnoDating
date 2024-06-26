plotPseudo = function(x) {
  if (class(r)!='resBactDating') error('plotPseudo needs a resBactDating object.')
  if (x$model!='poisson') error('plotPseudo only implements Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  ma=max(xs)*1.05
  rate=mean(x$record[(nrow(x$record)/2):nrow(x$record),'mu'])
  sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  par(mfrow=c(1,2))
  plot(c(0,ma),c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,ma),ylim=c(0,max(ys)*1.05))
  par(xpd=F)
  xss=seq(0,ma,ma/1000)
  plim=0.05
  
  lines(xss,qpois(  plim/2,xss*rate),lty='dashed')
  lines(xss,qpois(1-plim/2,xss*rate),lty='dashed')
  ll=dpois(round(ys),xs*rate,log=T)

  normed=(ll-min(ll))/(max(ll)-min(ll))
  cols=grDevices::rgb(1-normed,0,normed)
  base=seq(1,0,-0.1)
  legend("topleft",cex=0.5,legend=sprintf('%.3f',exp(base*(max(ll)-min(ll))+min(ll))),pch=19,col=grDevices::rgb(1-base,0,base))
  par(xpd=NA)
  points(xs,ys,pch=19,col=cols)
  par(xpd=F)
  
  plot(x$tree,show.tip.label = F,edge.color=cols)
  axisPhylo(1,backward = F)
}

plotPseudoQQ = function(x) {
  if (class(r)!='resBactDating') error('plotPseudo needs a resBactDating object.')
  if (x$model!='poisson') error('plotPseudo only implements Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=mean(x$record[(nrow(x$record)/2):nrow(x$record),'mu'])
  sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  p=ppois(round(ys),xs*rate)#uniform pseudo-residual
  p=runif(length(ys),ppois(round(ys-1),xs*rate),ppois(round(ys),xs*rate))
  n=qnorm(p)#normal pseudo-residual  
  mi=min(min(n),-3)
  ma=max(max(n),3)
  par(mfrow=c(2,2))
  hist(p,xlab='',main='Uniform pseudo-residuals',freq=F)
  lines(c(0,1),c(1,1))
  
  plot(n,xlab='',ylab='',main='Normal pseudo-residuals',ylim=c(mi,ma))
  lines(length(n)*c(-0.1,1.1),c(0,0))
  lines(length(n)*c(-0.1,1.1),rep(qnorm(0.005),2))
  lines(length(n)*c(-0.1,1.1),rep(qnorm(0.025),2))
  lines(length(n)*c(-0.1,1.1),rep(qnorm(0.975),2))
  lines(length(n)*c(-0.1,1.1),rep(qnorm(0.995),2))

  hist(n,xlab='',main='Normal pseudo-residuals',freq=F,xlim=c(mi,ma))
  xs=seq(mi,ma,(ma-mi)/1000)
  lines(xs,dnorm(xs))
  qqnorm(n,xlim=c(mi,ma),ylim=c(mi,ma))
  lines(c(mi,ma),c(mi,ma))
}