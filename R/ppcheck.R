#' Perform posterior predictive check
#'
#' @param x object of class resDating
#' @param nrep number of repeats to perform
#' @param showProgress Whether or not to show progress
#' @param showPlot Whether or not to show the plot
#' @export
#'
ppcheck=function(x,nrep=1000,showProgress=T,showPlot=F)
{
  phy=x$inputtree
  sta=computeStats(phy)
  if (!is.null(x$record)) inds=round(seq(max(1,floor(nrow(x$record)/2)),nrow(x$record),length.out=nrep))
  if (showProgress) pb <- utils::txtProgressBar(min=0,max=nrep,style = 3)
  staM=matrix(NA,nrow=nrep,ncol=length(sta))
  for (i in 1:nrep) {
    if (showProgress) utils::setTxtProgressBar(pb, i)
    if (is.null(x$record)) x2=x else x2=takeSample(x,inds[i])
    phy2=simobsphy(x2$tree,model=x2$model,mu=x2$rate,sigma=x2$relax)
    staM[i,]=computeStats(phy2)
  }

  if (showPlot){
    old.par=par(no.readonly = T)
    sub=ceiling(sqrt(length(sta)))
    par(mfrow=c(sub,sub))
    for (i in 1:length(sta)) {
      xl=c(min(staM[,i],sta[i]),max(staM[,i],sta[i]))
      hist(staM[,i],main='',xlab='',ylab='',xlim=xl)
      lines(c(sta[i],sta[i]),c(0,nrep),col='red')
      title(names(sta)[i])
    }
    par(old.par)
  }

  pvals=rep(NA,length(sta))
  for (i in 1:length(sta)) {
    pvals[i]=length(which(staM[,i]<sta[i]))/nrep
    if (pvals[i]>0.5) pvals[i]=1-pvals[i]#two-sided test
  }
  pvals=p.adjust(pvals,'fdr')
  return(min(pvals))
}

#Compute the summary statistics for a given phylogeny
computeStats=function(phy)
{
  phy=unroot(phy)
  sta=rep(NA,2)
  sta[1]=mean(phy$edge.length)#Mean of branch lengths
  sta[2]=var(phy$edge.length)#Variance of branch lengths
  sta[3]=max(phy$edge.length)#Max branch length
  sta[4]=sum(phy$edge.length[which(phy$edge[,2]>Ntip(phy))])/sta[1]#Stemminess
  names(sta)<-c('Mean bralen','Var bralen','Max bralen','Stemminess')
  return(sta)
}
