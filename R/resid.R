#' Calculate likelihood of branches
#'
#' @param x object of class resDating
#' @param log whether to return the log of the likelihoods
#'
calcLikBranches = function(x,log=FALSE) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=x$rate
  relax=x$relax

  if (x$model=='poisson') {
    probs=dpois(round(ys),xs*rate,log=log)
    return(probs)
  }

  if (x$model=='arc') {
    probs=dnbinom(round(ys),size=xs*rate/relax,prob=1/(1+relax),log=log)
    return(probs)
  }

  if (x$model=='strictgamma' || x$model=='carc') {
    probs=dgamma(ys,shape=xs*rate/(1+relax),scale=1+relax,log=log)
    return(probs)
  }

  stop(sprintf('Model %s is not yet implemented.',x$model))
}

#' Calculate normal residuals of branches
#'
#' @param x object of class resDating
#'
calcResiduals = function(x) {
  if (!inherits(x,'resDating')) stop('Not a resDating object.')
  xs=x$tree$edge.length
  ys=x$tree$subs
  rate=x$rate
  relax=x$relax

  if (x$model=='poisson') {
    probs=runif(length(ys),ppois(round(ys-1),xs*rate),ppois(round(ys),xs*rate))
    n=qnorm(probs)
    w=which(probs==1)
    if (length(w)>0) n[w]=qnorm(runif(length(w),ppois(round(ys[w]),xs[w]*rate,lower.tail = F),ppois(round(ys[w]-1),xs[w]*rate,lower.tail = F)),lower.tail = F)
    return(n)
  }

  if (x$model=='arc') {
    probs=runif(length(ys),pnbinom(round(ys-1),size=xs*rate/relax,prob=1/(1+relax)),pnbinom(round(ys),size=xs*rate/relax,prob=1/(1+relax)))
    n=qnorm(probs)
    w=which(probs==1)
    if (length(w)>0) n[w]=qnorm(runif(length(w),pnbinom(round(ys[w]),size=xs[w]*rate/relax,prob=1/(1+relax),lower.tail = F),pnbinom(round(ys[w]-1),size=xs[w]*rate/relax,prob=1/(1+relax),lower.tail = F)),lower.tail = F)
    return(n)
  }

  if (x$model=='strictgamma' || x$model=='carc') {
    probs=pgamma(ys,shape=xs*rate/(1+relax),scale=1+relax)
    n=qnorm(probs)
    w=which(probs==1)
    if (length(w)>0) n[w]=qnorm(pgamma(ys[w],shape=xs[w]*rate/(1+relax),scale=1+relax,lower.tail=F),lower.tail = F)
    return(n)
  }

  stop(sprintf('Model %s is not yet implemented.',x$model))
}


#' Plot likelihood of branches
#'
#' @param x object of class resDating
#' @param sub Plot only one of the four subplots
#' @param color Whether to use colors to show likelihoods
#' @param minProb Minimum likelihood to show as red
#' @param ... Passed on
#' @export
#'
plotLikBranches = function(x,sub=NA,color=T,minProb=NA,...) {
  ll=calcLikBranches(x,log=T)
  ll[is.infinite(ll)]=NA
  ll=pmin(ll,log(1))
  if (!is.na(minProb)) ll=pmax(ll,log(minProb))

  xs=x$tree$edge.length
  ys=x$tree$subs
  ma=max(xs)*1.05
  rate=x$rate
  relax=x$relax
  if (is.na(sub)) old.par=par(no.readonly = T)
  if (is.na(sub)) par(mfrow=c(1,2))
  if ((is.na(sub)) || sub==1) {
    def_args=list(x=c(0,ma),y=c(0,rate*ma),type='l',xlab='Branch duration',ylab='Substitutions',xaxs='i',yaxs='i',xlim=c(0,ma),ylim=c(0,max(ys)*1.05))
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub','color','minProb'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("plot",args)
    par(xpd=F)
    xss=seq(0,ma,ma/1000)
    plim=0.05

    if (x$model=='poisson') {
      lines(xss,qpois(  plim/2,xss*rate),lty='dashed')
      lines(xss,qpois(1-plim/2,xss*rate),lty='dashed')
    }

    if (x$model=='arc') {
      lines(xss,qnbinom(  plim/2,size=xss*rate/relax,prob=1/(1+relax)),lty='dashed')
      lines(xss,qnbinom(1-plim/2,size=xss*rate/relax,prob=1/(1+relax)),lty='dashed')
    }

    if (x$model=='strictgamma' || x$model=='carc') {
      lines(xss,qgamma(  plim/2,shape=xss*rate/(1+relax),scale=1+relax),lty='dashed')
      lines(xss,qgamma(1-plim/2,shape=xss*rate/(1+relax),scale=1+relax),lty='dashed')
    }
  }

  normed=(ll-min(ll,na.rm=T))/(max(ll,na.rm=T)-min(ll,na.rm=T))
  normed2=normed;normed2[is.na(normed2)]=1
  cols=ifelse(is.na(normed),grDevices::rgb(0,1,0),grDevices::rgb(1-normed2,0,normed2))
  base=seq(1,0,-0.1)
  if (color==F) cols=rep('black',length(cols))

  if ((is.na(sub)) || sub==1) {
    if (color) legend("topleft",cex=0.5,legend=sprintf('%.3f',exp(base*(max(ll,na.rm=T)-min(ll,na.rm=T))+min(ll,na.rm=T))),pch=19,col=grDevices::rgb(1-base,0,base))
    par(xpd=NA)
    points(xs,ys,pch=19,col=cols)
    par(xpd=F)
  }

  if (is.na(sub) || sub==2) {
    def_args=list(x=x$tree,show.tip.label = F,edge.color=cols)
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub','color','minProb'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("plot",args)
    axisPhylo(1,backward = F)
  }
  if (is.na(sub)) par(old.par)
}

#' Plot residuals
#'
#' @param x object of class resDating
#' @param sub Plot only one of the four subplots
#' @param ... Passed on
#' @export
#'
plotResid = function(x,sub=NA,...) {
  n=x$resid#normal residuals
  p=pnorm(n)#uniform residuals
  xs=x$tree$edge.length
  if (any(is.nan(n) | is.infinite(n))) {
    w=which(!is.nan(n) & !is.infinite(n))
    n=n[w]
    p=p[w]
    xs=xs[w]
    warning('Ignoring impossible branches.')
  }
  mi=min(min(n),-3)
  ma=max(max(n),3)
  if (is.na(sub)) old.par=par(no.readonly = T)
  if (is.na(sub)) par(mfrow=c(2,2))
  if (is.na(sub) || sub==1) {
    def_args=list(x=p,xlab='',main='Uniform residuals',freq=F)
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("hist",args)
    abline(1,0,lty=3)
  }
  o=order(xs)

  if (is.na(sub) || sub==2) {
    def_args=list(x=n[o],xlab='Branches in increasing order of duration',ylab='',main='Normal residuals',ylim=c(mi,ma))
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("plot",args)
    abline(h=qnorm(c(0.005,0.025,0.5,0.975,0.995)),lty=3)
  }

  if (is.na(sub) || sub==3) {
    def_args=list(x=n,xlab='',main='Normal residuals',freq=F,xlim=c(mi,ma))
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("hist",args)
    xs=seq(mi,ma,(ma-mi)/1000)
    par(xpd=NA)
    lines(xs,dnorm(xs),lty=3)
    par(xpd=F)
  }

  if (is.na(sub) || sub==4) {
    P <- ppoints(length(n))
    z=qnorm(P)
    def_args=list(x=z,y=sort(n),xlab='Theoretical quantiles',ylab='Sample quantiles', main = "Normal Q-Q Plot",xlim=c(mi,ma),ylim=c(mi,ma))
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("plot",args)
    P1000 = ppoints(1000)
    z=qnorm(P1000)
    grid(lty=1)
    abline(0,1,lty=3)
    SE <- (1/dnorm(z))*sqrt(P1000*(1 - P1000)/length(n))
    upper <- z+SE*1.959964
    lower <- z-SE*1.959964
    lines(z, upper, lty=3)
    lines(z, lower, lty=3)
  }
  if (is.na(sub)) par(old.par)
}

#' Test on residuals
#'
#' @param x object of class resDating
#' @param test which test to use: 1 for Anderson-Darling simple test (default), 2 for Kolmogorov-Smirnov simple test, 3 for Shapiro-Wilk composite test
#' @export
#'
testResid=function(x,test=1) {
  n=x$resid#normal residual
  p=pnorm(n)#uniform residual
  if (any(is.nan(n) | is.infinite(n))) {
    w=which(!is.nan(n) & !is.infinite(n))
    n=n[w]
    p=p[w]
    warning('Ignoring impossible branches.')
  }
  if (test==1) r=DescTools::AndersonDarlingTest(n,null='pnorm',mean=0,sd=1)
  #the above is same as DescTools::AndersonDarlingTest(p,null='punif',min=0,max=1)
  #and also same as ADGofTest::ad.test(n,pnorm)
  if (test==2) r=ks.test(n,pnorm,mean=0,sd=1)
  #the above is same as ks.test(p,punif,min=0,max=1)
  if (test==3) r=shapiro.test(n)
  return(r)
}

#' Resample a dated tree using bactdate
#'
#' @param x object of class resDating
#' @param showProgress Whether or not to show progress
#' @param showTraces whether or not to show the MCMC traces
#' @param nbIts number of iterations in bactdate resampling
#' @export
#'
resample=function(x,showProgress=T,showTraces=F,nbIts=1e4)
{
  rtree=x$tree#reorder(x$tree)
  phy=x$inputtree
  dates=unname(x$rootdate+dist.nodes(x$tree)[1:Ntip(x$tree),1+Ntip(x$tree)])
  rtree$root.time=NULL
  k=sum(phy$edge.length)/sum(x$tree$edge.length)
  rtree$edge.length=rtree$edge.length*k
  rtree$edge.length=ifelse(rtree$edge.length==0,0.01,rtree$edge.length)
  r2=bactdate(rtree,dates,initMu=k,minbralen=0,model='strictgamma',showProgress=showProgress,initAlpha=estimAlpha(x$tree),updateRoot = 'branch',nbIts=nbIts)
  attributes(r2$tree)$order<-NULL
  attributes(r2$inputtree)$order<-NULL
  if (showTraces) {tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')}
  r4=resDating(r2$tree,phy,algo=x$algo,model=x$model,rate=x$rate,relax=x$relax)
  r4$record=r2$record
  x3=resDating(takeSample(r4)$tree,r4$inputtree,algo=r4$algo,model=r4$model,rate=r4$rate,relax=r4$relax,rootdate=r4$rootdate)
  s=Ntip(rtree)+Nnode(rtree)
  i2=1:nrow(r4$tree$edge)
  for (w in 1:nrow(r4$record)) r4$record[w,s*2+r4$tree$edge[i2,2]]=x3$tree$subs
  r4$record[,'mu']=x$rate
  r4$record[,'sigma']=x$relax
  return(r4)
}

#' Compute posterior distribution of p-values
#'
#' @param x object of class resDating
#' @param nrep number of repeats to perform
#' @export
#'
postdistpvals=function(x,nrep=1000)
{
  resampling=0;showProgress=T;shoeTraces=F;nbIts=1e4#old arguments
  ps=rep(NA,nrep)

  if (resampling==0) {
    #Using existing posterior sample
    if (is.null(x$record)) stop('No posterior sample was found.')
    inds=round(seq(max(1,floor(nrow(x$record)/2)),nrow(x$record),length.out=nrep))
    for (i in 1:nrep) {
      x2=takeSample(x,inds[i])
      ps[i]=testResid(x2)$p.value
    }
  }

  if (resampling==2) {
    #Using resampling via bactdate
    rtree=x$tree#reorder(x$tree)
    phy=x$inputtree
    dates=unname(x$rootdate+dist.nodes(x$tree)[1:Ntip(x$tree),1+Ntip(x$tree)])
    rtree$root.time=NULL
    k=sum(phy$edge.length)/sum(x$tree$edge.length)#or use x$rate
    rtree$edge.length=rtree$edge.length*k
    rtree$edge.length=ifelse(rtree$edge.length==0,0.01,rtree$edge.length)
    r2=bactdate(rtree,dates,initMu=k,minbralen=0,model='strictgamma',showProgress=showProgress,initAlpha=estimAlpha(x$tree),updateRoot = 'branch',nbIts=nbIts)
    attributes(r2$tree)$order<-NULL
    attributes(r2$inputtree)$order<-NULL
    if (showTraces) {tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')}
    r4=resDating(r2$tree,phy,algo=x$algo,model=x$model,rate=x$rate,relax=x$relax)
    r4$record=r2$record
    x3=resDating(takeSample(r4)$tree,r4$inputtree,algo=r4$algo,model=r4$model,rate=r4$rate,relax=r4$relax,rootdate=r4$rootdate)
    inds=round(seq(max(1,floor(nrow(r4$record)/2)),nrow(r4$record),length.out=nrep))
    for (i in 1:nrep) {
      x2=takeSample(r4,inds[i])
      x2$tree$subs=x3$tree$subs
      x2$rate=x$rate
      x2$relax=x$relax
      x2$resid=calcResiduals(x2)
      ps[i]=testResid(x2)$p.value
    }
  }

  if (resampling==3) {
    #Using resampling via bactdate, slower but clearer version
    rtree=x$tree#reorder(x$tree)
    phy=x$inputtree
    dates=unname(x$rootdate+dist.nodes(x$tree)[1:Ntip(x$tree),1+Ntip(x$tree)])
    rtree$root.time=NULL
    k=x$rate#sum(phy$edge.length)/sum(rtree$edge.length)
    rtree$edge.length=rtree$edge.length*k
    r2=runDating(rtree,dates,rate=k,minbralen=0,algo='BactDating',model='strictgamma',showProgress=showProgress,updateRoot='branch',initAlpha=estimAlpha(x$tree),nbIts=nbIts)
    if (showTraces) {tmp=r2;class(tmp)='resBactDating';plot(tmp,'trace')}
    inds=round(seq(max(1,floor(nrow(r2$record)/2)),nrow(r2$record),length.out=nrep))
    if (showProgress) print('Computing p-values...')
    if (showProgress) pb <- utils::txtProgressBar(min=0,max=nrep,style = 3)
    for (i in 1:nrep) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      samtree=takeSample(r2,inds[i])$tree
      x2=resDating(samtree,phy,algo=x$algo,model=x$model,rate=x$rate,relax=x$relax)
      ps[i]=testResid(x2)$p.value
    }
    if (showProgress) close(pb)
  }

  if (resampling==1) {
    #Generating approximate posterior sample
    dates=unname(dist.nodes(x$tree)[Ntip(x$tree)+1,1:Ntip(x$tree)])
    mu=x$rate
    l2=unname(x$tree$edge.length)
    n=length(l2)
    store=matrix(NA,1000,n)
    for (i in 1:nrow(store)) {
      t=simcoaltree(dates,alpha=sampleAlpha(x$tree))$edge.length
      s2=rpois(n,t*mu)
      ind=order(s2)
      t=t[ind]
      s2=s2[ind]
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
    }
    x2=x
    for (i in 1:nrep) {
      x2$tree$edge.length=rgamma(n,shape=k+l2*mu,scale=theta/(1+theta*mu))
      x2$resid=calcResiduals(x2)
      ps[i]=testResid(x2)$p.value
    }
  }

  ret=list()
  ret$ps=ps
  ret$input=x
  ret$last=x2
  class(ret)<-'resPDPV'
  return(ret)
}

#' Print function for resPDPV objects
#' @param x output from postdistpval
#' @param ... Passed on to cat
#' @return Print out details of dating results
#' @export
print.resPDPV <- function(x, ...)
{
  stopifnot(inherits(x, "resPDPV"))
  cat(sprintf('Result from postdistpval with median p-value %.2f\n',median(x$ps)),...)
  invisible(x)
}

#' Plotting methods
#' @param x Output from postdistpval
#' @param type Type of plot to do. Currently either 'hist' or 'lastLikBranches' or 'lastResid'.
#' @param ... Additional parameters are passed on
#' @return Plot of results
#' @export
plot.resPDPV = function(x, type='hist',...) {
  stopifnot(inherits(x, "resPDPV"))

  if (type=='lastLikBranches') {
    plotLikBranches(x$last,...)
  }

  if (type=='lastResid') {
    plotResid(x$last,...)
    title(sprintf('p=%s',format(testResid(x$last)$p.value)))
    }

  if (type=='hist') {
    h=hist(x$ps,breaks=seq(0,1,length.out=21),xlab='',ylab='',main='Posterior distribution of p-values')
    y=length(which(x$ps<0.05))
    par(xpd=NA)
    text(0.025,y,sprintf('%.1f%%',y*100/length(x$ps)),pos=3)
    med=median(x$ps)
    lines(c(med,med),c(0,max(h$counts)),lty=2)
    text(med,max(h$counts),format(med,digits=3),pos=4)
    par(xpd=F)
  }
}
