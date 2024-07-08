#' Calculate probability of branches
#'
#' @param x object of class resDating
#' @param cumul whether to return the cumulative probability
#' @param log whether to return the log of the probability
#'
calcProbBranches = function(x,cumul=FALSE,log=FALSE) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  if (x$model!='poisson') stop('Only Poisson model at the moment.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=x$rate
  #sigma=mean(x$record[(nrow(x$record)/2):nrow(x$record),'sigma'])
  if (cumul==FALSE) {
    probs=dpois(round(ys),xs*rate,log=log)
  } else {
    #probs=ppois(round(ys),xs*rate,log.p=log)
    oldseed <- .Random.seed
    set.seed(1998)
    probs=runif(length(ys),ppois(round(ys-1),xs*rate),ppois(round(ys),xs*rate))
    .Random.seed <- oldseed
    if (log) probs=log(probs)
  }
  return(probs)
}

#' Plot probability of branches
#'
#' @param x object of class resDating
#' @export
#'
plotProbBranches = function(x) {
  ll=calcProbBranches(x,log=T)
  ll[is.infinite(ll)]=NA

  xs=x$tree$edge.length
  ys=x$tree$subs
  ma=max(xs)*1.05
  rate=x$rate
  old.par=par(no.readonly = T)
  par(mfrow=c(1,2))
  plot(c(0,ma),c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,ma),ylim=c(0,max(ys)*1.05))
  par(xpd=F)
  xss=seq(0,ma,ma/1000)
  plim=0.05

  lines(xss,qpois(  plim/2,xss*rate),lty='dashed')
  lines(xss,qpois(1-plim/2,xss*rate),lty='dashed')

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
  par(old.par)
}

#' Plot pseudo-residuals
#'
#' @param x object of class resDating
#' @export
#'
plotResid = function(x) {
  p=calcProbBranches(x,cumul=T)#uniform pseudo-residual
  if (any(p==1)) {
    p=p[which(p!=1)]
    warning('Ignoring branches with zero probability.')
  }
  n=qnorm(p)#normal pseudo-residual
  mi=min(min(n),-3)
  ma=max(max(n),3)
  old.par=par(no.readonly = T)
  par(mfrow=c(2,2))
  hist(p,xlab='',main='Uniform pseudo-residuals',freq=F)
  abline(1,0,lty=3)

  o=order(x$tree$edge.length)
  plot(n[o],xlab='Branches in increasing order of duration',ylab='',main='Normal pseudo-residuals',ylim=c(mi,ma))
  abline(h=qnorm(c(0.005,0.025,0.5,0.975,0.995)),lty=3)

  hist(n,xlab='',main='Normal pseudo-residuals',freq=F,xlim=c(mi,ma))
  xs=seq(mi,ma,(ma-mi)/1000)
  par(xpd=NA)
  lines(xs,dnorm(xs),lty=3)
  par(xpd=F)

  P <- ppoints(length(n))
  z=qnorm(P)
  plot(z, sort(n),xlab='Theoretical quantiles',ylab='Sample quantiles', main = "Normal Q-Q Plot",xlim=c(mi,ma),ylim=c(mi,ma))
  P1000 = ppoints(1000)
  z=qnorm(P1000)
  grid(lty=1)
  abline(0,1,lty=3)
  SE <- (1/dnorm(z))*sqrt(P1000*(1 - P1000)/length(n))
  upper <- z+SE*1.959964
  lower <- z-SE*1.959964
  lines(z, upper, lty=3)
  lines(z, lower, lty=3)

  par(old.par)
}

#' Test on pseudo-residuals
#'
#' @param x object of class resDating
#' @export
#'
testResid=function(x) {
  p=calcProbBranches(x,cumul=T)#uniform pseudo-residual
  if (any(p==1)) {
    p=p[which(p!=1)]
    warning('Ignoring branches with zero probability.')
  }
  n=qnorm(p)#normal pseudo-residual
  r1=shapiro.test(n)#tests if Normal but not if Normal(0,1)
  r2=ks.test(p,punif,min=0,max=1)#this is same as ks.test(n,pnorm,mean=0,sd=1)
  if (r1$p.value<r2$p.value) return(r1) else return(r2)
}
