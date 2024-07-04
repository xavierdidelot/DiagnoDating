#' Run root-to-tip analysis
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to BactDating::roottotip
#'
#' @return resDating object containing results of root-to-tip analysis
#' @export
#'
runRTT=function(tree,dates,...) {
  suppressWarnings(r<-BactDating::roottotip(tree,dates,...))
  r$algo='RTT'
  r$rootdate=r$ori
  class(r)<-'resDating'
  return(r)
}

#' Date a tree using BactDating
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to BactDating::bactdate
#'
#' @return resDating object containing results of BactDating analysis
#' @export
#'
runBactDating=function(tree,dates,...) {
  r=BactDating::bactdate(tree,dates,...)
  r$algo='BactDating'
  v=r$record[,'mu']
  v=v[(1+length(v)/2):length(v)]
  r$rate=mean(v)
  class(r)<-'resDating'
  return(r)
}

#' Date a tree using treedater
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to treedater::dater
#'
#' @return resDating object containing results of treedater analysis
#' @export
#'
runTreeDater=function(tree,dates,...) {
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l,...)))
  res=list()
  res$algo='treedater'
  res$inputtree=tree
  res$model='poisson'
  res$rate=rtd$mean.rate*l
  res$rootdate=rtd$timeOfMRCA
  nrowtab=Ntip(tree)+Nnode(tree)
  record = matrix(NA, 10, nrowtab*3 + 6)
  colnames(record)<-c(rep(NA,nrowtab*3),'likelihood','mu','sigma','alpha','prior','root')
  record[,'likelihood']=rtd$loglik
  record[,'mu']=res$rate
  record[,Ntip(tree)+1]=rtd$timeOfMRCA
  res$record=record
  res$tree=list(edge=rtd$edge,Nnode=rtd$Nnode,tip.label=rtd$tip.label,edge.length=rtd$edge.length,root.time=rtd$timeOfMRCA)
  class(res$tree)<-'phylo'
  tree=reorderEdges(tree,res$tree)
  res$tree$subs=tree$edge.length
  res$algo='treeedater'
  class(res)<-'resDating'
  return(res)
}

#' Date a tree using LSD
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to Rlsd2::lsd2
#'
#' @return resDating object containing results of LSD analysis
#' @export
#'
runLSD=function(tree,dates,...) {
  tag=round(runif(1,1,1e8))
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.csv',tag),quote = F,col.names=length(sts))
  o=capture.output(result <- Rlsd2::lsd2(inputTree=sprintf('/tmp/tree%d.nwk',tag), inputDate=sprintf('/tmp/dates%d.csv',tag),outFile = sprintf('/tmp/result%d',tag), seqLen=l,nullblen=-1,...))
  res=list()
  res$inputtree=tree
  res$model='poisson'
  res$rate=result$rate*l
  res$rootdate=result$tMRCA
  nrowtab=Ntip(tree)+Nnode(tree)
  record = matrix(NA, 10, nrowtab*3 + 6)
  colnames(record)<-c(rep(NA,nrowtab*3),'likelihood','mu','sigma','alpha','prior','root')
  record[,'mu']=result$rate*l
  record[,Ntip(tree)+1]=result$tMRCA
  res$record=record
  res$tree=read.nexus(sprintf('/tmp/result%d.date.nexus',tag))
#  res$tree=multi2di(res$tree)
  tree=reorderEdges(tree,res$tree)
  res$tree$subs=tree$edge.length
  res$tree$root.time=res$rootdate
  res$algo='LSD'
  class(res)<-'resDating'
  return(res)
}

#' Print function for resDating objects
#' @param x output from bactdate
#' @param ... Passed on to cat
#' @return Print out details of dating results
#' @export
print.resDating <- function(x, ...)
{
  stopifnot(inherits(x, "resDating"))
  cat( 'Result of analysis using',x$algo,'- clock rate is',x$rate, 'and root date is',x$rootdate,...)
  invisible(x)
}

#' Plotting methods
#' @param x Output from dating
#' @param ... Additional parameters are passed on
#' @return Plot of results
#' @export
plot.resDating = function(x, ...) {
  if (x$algo=='RTT') stop('Plotting not available.')
  class(x)<-'resBactDating'
  plot(x,...)
}
