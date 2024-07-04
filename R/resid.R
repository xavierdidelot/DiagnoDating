#' Plot probability of branches
#'
#' @param x object of class resDating
#' @export
#'
plotProbBranches = function(x) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  if (x$model!='poisson') stop('Only Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  ma=max(xs)*1.05
  rate=x$rate
  #sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  par(mfrow=c(1,2))
  plot(c(0,ma),c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,ma),ylim=c(0,max(ys)*1.05))
  par(xpd=F)
  xss=seq(0,ma,ma/1000)
  plim=0.05

  lines(xss,qpois(  plim/2,xss*rate),lty='dashed')
  lines(xss,qpois(1-plim/2,xss*rate),lty='dashed')
  ll=dpois(round(ys),xs*rate,log=T)
  ll[is.infinite(ll)]=NA

  normed=(ll-min(ll,na.rm=T))/(max(ll,na.rm=T)-min(ll,na.rm=T))
  normed2=normed;normed2[is.na(normed2)]=1
  cols=ifelse(is.na(normed),grDevices::rgb(0,1,0),grDevices::rgb(1-normed2,0,normed2))
  base=seq(1,0,-0.1)
  legend("topleft",cex=0.5,legend=sprintf('%.3f',exp(base*(max(ll,na.rm=T)-min(ll,na.rm=T))+min(ll,na.rm=T))),pch=19,col=grDevices::rgb(1-base,0,base))
  par(xpd=NA)
  points(xs,ys,pch=19,col=cols)
  par(xpd=F)

  plot(x$tree,show.tip.label = F,edge.color=cols)
  axisPhylo(1,backward = F)
}

#' Plot pseudo-residuals
#'
#' @param x object of class resDating
#' @export
#'
plotResid = function(x) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  if (x$model!='poisson') stop('Only Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=x$rate
  #sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  p=ppois(round(ys),xs*rate)#uniform pseudo-residual
  if (any(p==1)) {
    p=p[which(p!=1)]
    warning('Removing p=1')
  }
  #p=(ppois(round(ys-1),xs*rate)+ppois(round(ys),xs*rate))/2
  #p=runif(length(ys),ppois(round(ys-1),xs*rate),ppois(round(ys),xs*rate))
  n=qnorm(p)#normal pseudo-residual
  n=(n-mean(n))/sd(n);p=pnorm(n)#weird normalize
  mi=min(min(n),-3)
  ma=max(max(n),3)
  par(mfrow=c(2,2))
  hist(p,xlab='',main='Uniform pseudo-residuals',freq=F)
  lines(c(0,1),c(1,1))

  plot(n,xlab='',ylab='',main='Normal pseudo-residuals',ylim=c(mi,ma))
  abline(h=qnorm(c(0.005,0.025,0.5,0.975,0.995)))

  hist(n,xlab='',main='Normal pseudo-residuals',freq=F,xlim=c(mi,ma))
  xs=seq(mi,ma,(ma-mi)/1000)
  lines(xs,dnorm(xs))

  qqnorm(n,xlim=c(mi,ma),ylim=c(mi,ma))
  abline(0,1)
}

#' Test on pseudo-residuals
#'
#' @param x object of class resDating
#' @export
#'
testResid=function(x) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  if (x$model!='poisson') stop('Only Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=x$rate
  #sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  p=ppois(round(ys),xs*rate)#uniform pseudo-residual
  if (any(p==1)) {
    p=p[which(p!=1)]
    warning('Removing p=1')
  }
  n=qnorm(p)#normal pseudo-residual
  r=shapiro.test(n)
  return(r)
}
