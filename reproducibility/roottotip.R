roottotip = function(tree,date,rate=NA,permTest=10000,showFig=T,colored=T,showPredInt='gamma',showText=T,showTree=T)
{
  if (!is.rooted(tree)) warning('Warning: roottotip was called on an unrooted input tree. Consider using initRoot first.\n')
  if (sum(tree$edge.length)<5) warning('Warning: input tree has small branch lengths. Make sure branch lengths are in number of substitutions (NOT per site).\n')
  #Rerranging of dates, if needed
  if (!is.null(names(date))) date=findDates(tree,date)

  if (var(date,na.rm=T)==0 && is.na(rate)) {warning('Warning: All dates are identical.\n');return(list(rate=NA,ori=NA,pvalue=NA))}
  n=length(date)
  ys=leafDates(tree)
  if (is.na(rate)) {
    res=lm(ys~date)
  }
  else {
    res=lm(I(ys-rate*date)~1)
    res$coefficients=c(res$coefficients,rate)
  }
  ori=-coef(res)[1]/coef(res)[2]
  rate=coef(res)[2]
  r2=summary(res)$r.squared
  correl=cor(date,ys,use='complete.obs')
  #pvalue=summary(res)$coefficients[,4][2]
  #print(c(r2,correl^2))#Equal

  pvalue=0
  for (i in 1:permTest) {
    date2=sample(date,n,replace=F)
    correl2=cor(date2,ys,use='complete.obs')
    if (correl2>=correl) pvalue=pvalue+1/permTest
  }

  if (rate<0) {warning('The linear regression suggests a negative rate.')}
  if (showFig==F) return(list(rate=rate,ori=ori,pvalue=pvalue))
  old.par=par(no.readonly = T)
  par(xpd=NA,oma = c(0, 0, 2, 0))
  if (colored) {
    normed=(date-min(date,na.rm=T))/(max(date,na.rm=T)-min(date,na.rm=T))
    cols=grDevices::rgb(ifelse(is.na(normed),0,normed),ifelse(is.na(normed),0.5,0),1-ifelse(is.na(normed),1,normed),0.5)
  } else cols='black'
  if (showTree) {
    #par(mfrow=c(1,2))
    plot(tree,show.tip.label = F)
    if (colored) tiplabels(col=cols,pch=19)
    axisPhylo(1,backward = F)
    fig_label('A',cex=2)
  }
  plot(date,ys,col=cols,xlab=ifelse(showText,'Sampling date',''),ylab=ifelse(showText,'Root-to-tip distance',''),xaxs='i',yaxs='i',pch=19,ylim=c(0,max(ys)),xlim=c(ifelse(rate>0,ori,min(date,na.rm = T)),max(date,na.rm = T)))
  #text(date,ys,labels=1:length(date))
  par(xpd=F)
  abline(res,lwd=2)
  if (rate<0) {par(old.par);return(list(rate=rate,ori=ori,pvalue=pvalue))}
  xs=seq(ori,max(date,na.rm = T),0.1)
  plim=0.05
  if (showPredInt=='poisson') {
    lines(xs,qpois(  plim/2,(xs-ori)*rate),lty='dashed')
    lines(xs,qpois(1-plim/2,(xs-ori)*rate),lty='dashed')
  }
  if (showPredInt=='gamma') {
    lines(xs,qgamma(  plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
    lines(xs,qgamma(1-plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
  }
  #if (showText) {
  #  if (pvalue==0) mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p<%.2e',rate,ori,r2,1/permTest), outer = TRUE, cex = 1.5)
  #  else           mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p=%.2e',rate,ori,r2,pvalue), outer = TRUE, cex = 1.5)
  #}
  #par(old.par)
  return(list(rate=rate,ori=ori,pvalue=pvalue))
}
