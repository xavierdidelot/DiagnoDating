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