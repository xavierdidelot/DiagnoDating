#' Calculate probability of branches
#'
#' @param x object of class resDating
#' @param log whether to return the log of the probability
#'
calcProbBranches = function(x,log=FALSE) {
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

#' Calculate residuals of branches
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
    return(probs)
  }

  if (x$model=='arc') {
    probs=runif(length(ys),pnbinom(round(ys-1),size=xs*rate/relax,prob=1/(1+relax)),pnbinom(round(ys),size=xs*rate/relax,prob=1/(1+relax)))
    return(probs)
  }

  if (x$model=='strictgamma' || x$model=='carc') {
    probs=pgamma(ys,shape=xs*rate/(1+relax),scale=1+relax)
    return(probs)
  }

  stop(sprintf('Model %s is not yet implemented.',x$model))
}


#' Plot probability of branches
#'
#' @param x object of class resDating
#' @param sub Plot only one of the four subplots
#' @param color Whether to use colors to show probabilities
#' @param minProb Minimum probability to show as red
#' @param ... Passed on
#' @export
#'
plotProbBranches = function(x,sub=NA,color=T,minProb=NA,...) {
  ll=calcProbBranches(x,log=T)
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

#' Plot pseudo-residuals
#'
#' @param x object of class resDating
#' @param sub Plot only one of the four subplots
#' @param ... Passed on
#' @export
#'
plotResid = function(x,sub=NA,...) {
  p=x$resid#uniform pseudo-residual
  xs=x$tree$edge.length
  if (any(is.nan(p) | p==0 | p==1)) {
    w=which(!is.nan(p) & p!=0 & p!=1)
    p=p[w]
    xs=xs[w]
    warning('Ignoring impossible branches.')
  }
  n=qnorm(p)#normal pseudo-residual
  mi=min(min(n),-3)
  ma=max(max(n),3)
  if (is.na(sub)) old.par=par(no.readonly = T)
  if (is.na(sub)) par(mfrow=c(2,2))
  if (is.na(sub) || sub==1) {
    def_args=list(x=p,xlab='',main='Uniform pseudo-residuals',freq=F)
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("hist",args)
    abline(1,0,lty=3)
  }
  o=order(xs)

  if (is.na(sub) || sub==2) {
    def_args=list(x=n[o],xlab='Branches in increasing order of duration',ylab='',main='Normal pseudo-residuals',ylim=c(mi,ma))
    cl=as.list(match.call())[-1]
    cl=cl[setdiff(names(cl),c('x','sub'))]
    args=c(cl,def_args[!names(def_args) %in% names(cl)])
    do.call("plot",args)
    abline(h=qnorm(c(0.005,0.025,0.5,0.975,0.995)),lty=3)
  }

  if (is.na(sub) || sub==3) {
    def_args=list(x=n,xlab='',main='Normal pseudo-residuals',freq=F,xlim=c(mi,ma))
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

#' Test on pseudo-residuals
#'
#' @param x object of class resDating
#' @param test which test to use: 1 for Anderson-Darling simple test (default), 2 for Kolmogorov-Smirnov simple test, 3 for Shapiro-Wilk composite test
#' @export
#'
testResid=function(x,test=1) {
  p=x$resid#uniform pseudo-residual
  if (any(is.nan(p) | p==0 | p==1)) {
    p=p[which(!is.nan(p) & p!=0 & p!=1)]
    warning('Ignoring impossible branches.')
  }
  n=qnorm(p)#normal pseudo-residual
  if (test==1) r=DescTools::AndersonDarlingTest(n,null='pnorm',mean=0,sd=1)
  #the above is same as DescTools::AndersonDarlingTest(p,null='punif',min=0,max=1)
  #and also same as ADGofTest::ad.test(n,pnorm)
  if (test==2) r=ks.test(n,pnorm,mean=0,sd=1)
  #the above is same as ks.test(p,punif,min=0,max=1)
  if (test==3) r=shapiro.test(n)
  return(r)
}

#' Check whether a dating analysis is valid
#'
#' @param x object of class resDating
#' @param nrep number of repeats to perform
#' @param resampling method used for resampling: 0 (default) for using given sample, 1 to use the approximate posterior resampling method, 2 to use BactDating resampling
#' @param nstore size of the store to use in approximate posterior resampling method
#' @param showProgress Whether or not to show progress
#' @param showPlot whether or not to show the histogram of p-values
#' @param showLast whether or not to show the residuals for the last sample
#' @param showTraces whether or not to show the MCMC traces
#' @export
#'
validate=function(x,nrep=1000,resampling=0,nstore=1000,showProgress=T,showPlot=T,showLast=F,showTraces=F)
{

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
    rtree=reorder(x$tree)
    phy=x$inputtree
    dates=unname(x$rootdate+dist.nodes(x$tree)[1:Ntip(x$tree),1+Ntip(x$tree)])
    rtree$root.time=NULL
    k=x$rate#sum(phy$edge.length)/sum(rtree$edge.length)
    rtree$edge.length=rtree$edge.length*k
    r2=runDating(rtree,dates,rate=k,minbralen=1e-10,algo='BactDating',model='strictgamma',updateRoot='branch',showProgress=showProgress,initAlpha=estimAlpha(x$tree),updateAlpha=F)
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
    rtree=reorder(x$tree)
    phy=x$inputtree
    dates=unname(x$rootdate+dist.nodes(x$tree)[1:Ntip(x$tree),1+Ntip(x$tree)])
    rtree$root.time=NULL
    k=x$rate#sum(phy$edge.length)/sum(rtree$edge.length)
    rtree$edge.length=rtree$edge.length*k
    r2=runDating(rtree,dates,rate=k,minbralen=1e-10,algo='BactDating',model='strictgamma',showProgress=showProgress,updateRoot='branch',initAlpha=estimAlpha(x$tree),updateAlpha=F)
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
  }

  if (resampling==1) {
    #Generating approximate posterior sample
    dates=unname(dist.nodes(x$tree)[Ntip(x$tree)+1,1:Ntip(x$tree)])
    mu=x$rate
    l2=unname(x$tree$edge.length)
    n=length(l2)
    store=matrix(NA,nstore,n)
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

  #Show the probability of branches and distribution of residuals in the last iteration if requested
  if (showLast) {plotProbBranches(x2);plotResid(x2);title(testResid(x2)$p.value)}

  #Plot histogram if requested
  if (showPlot) {
    h=hist(ps,breaks=seq(0,1,length.out=21),xlab='',ylab='',main='Posterior distribution of p-values')
    y=length(which(ps<0.05))
    par(xpd=NA)
    text(0.025,y,sprintf('%.1f%%',y*100/nrep),pos=3)
    med=median(ps)
    lines(c(med,med),c(0,max(h$counts)),lty=2)
    text(med,max(h$counts),format(med,digits=3),pos=4)
  }
  invisible(ps)
}
