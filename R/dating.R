#' Run a dating analysis
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param algo Algorithm to use, can be one of: LSD, node.dating, BactDating, treedater, TreeTime
#' @param ... Passed on to dating algorithm
#'
#' @return resDating object containing results of dating analysis
#' @export
#'
runDating=function(tree,dates,algo='BactDating',...) {
  if (algo=='1' || algo=='LSD') r=runLSD(tree,dates,...)
  if (algo=='2' || algo=='node.dating') r=runNodeDating(tree,dates,...)
  if (algo=='3' || algo=='BactDating') r=runBactDating(tree,dates,...)
  if (algo=='4' || algo=='treedater') r=runTreeDater(tree,dates,...)
  if (algo=='5' || algo=='TreeTime') r=runTreeTime(tree,dates,...)
  return(r)
}

#' Date a tree using BactDating
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to BactDating::bactdate
#'
#' @return resDating object containing results of BactDating analysis
#'
runBactDating=function(tree,dates,...) {
  def_args=list(tree=tree,date=dates,model='poisson',updateRoot=F)
  cl=as.list(match.call())[-(1:3)]
  args=c(cl,def_args[!names(def_args) %in% names(cl)])
  r=do.call("bactdate",args)
  r$algo='BactDating'
  r$model=args$model
  v=r$record[,'mu']
  v=v[(1+length(v)/2):length(v)]
  r$rate=mean(v)
  v=r$record[,'sigma']
  v=v[(1+length(v)/2):length(v)]
  r$relax=mean(v)
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
#'
runTreeDater=function(tree,dates,...) {
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tree$tip.label
  o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l,...)))
  res=list()
  res$algo='treedater'
  res$inputtree=tree
  res$model='poisson'
  res$rate=rtd$mean.rate*l
  res$relax=0
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
#'
runLSD=function(tree,dates,...) {
  tag=round(runif(1,1,1e8))
  tre=tree
  l=round(max(tree$edge.length)*1000)
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tre$tip.label
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.csv',tag),quote = F,col.names=length(sts))
  o=capture.output(result <- Rlsd2::lsd2(inputTree=sprintf('/tmp/tree%d.nwk',tag), inputDate=sprintf('/tmp/dates%d.csv',tag),outFile = sprintf('/tmp/result%d',tag), seqLen=l,nullblen=-1,...))
  res=list()
  res$inputtree=tree
  res$model='poisson'
  res$rate=result$rate*l
  res$relax=0
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

#' Date a tree using node.dating
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on to ape::node.dating
#'
#' @return resDating object containing results of node.dating analysis
#'
runNodeDating=function(tree,dates,...) {
  try(suppressWarnings(mu<-ape::estimate.mu(tree,dates)),silent=T)
  if (mu<=0) mu=0.01
  suppressWarnings(d<-ape::estimate.dates(tree,dates,mu=mu))
  res=list()
  res$algo='node.dating'
  res$model='poisson'
  res$rate=mu
  res$relax=0
  res$rootdate=min(d)
  dt=tree
  dt$subs=dt$edge.length
  dt$root.time=res$rootdate
  for (i in 1:nrow(dt$edge)) dt$edge.length[i]=d[dt$edge[i,2]]-d[dt$edge[i,1]]
  res$inputtree=tree
  res$tree=dt
  class(res)<-'resDating'
  return(res)
}

#' Date a tree using TreeTime
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param ... Passed on
#'
#' @return resDating object containing results of TreeTime analysis
#'
runTreeTime=function(tree,dates,...) {
  tag=round(runif(1,1,1e8))
  tre=tree
  l=round(max(tree$edge.length)*1000)
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tree$tip.label
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.tsv',tag),quote = F,col.names='strain\tdate',sep='\t')
  system(sprintf("treetime --tree /tmp/tree%d.nwk --dates /tmp/dates%d.tsv --sequence-length %d --keep-root --outdir /tmp/%d",tag,tag,l,tag))
  res=list()
  res$inputtree=tree
  res$model='poisson'
  res$rate=read.table(sprintf('/tmp/%d/molecular_clock.txt',tag))[1,1]*l
  res$relax=0
  res$tree=read.nexus(sprintf('/tmp/%d/timetree.nexus',tag))
  res$rootdate=max(dates)-max(dist.nodes(res$tree)[Ntip(res$tree)+1,1:Ntip(res$tree)])
  nrowtab=Ntip(tree)+Nnode(tree)
  record = matrix(NA, 10, nrowtab*3 + 6)
  colnames(record)<-c(rep(NA,nrowtab*3),'likelihood','mu','sigma','alpha','prior','root')
  record[,'mu']=res$rate
  record[,Ntip(tree)+1]=res$rootdate
  res$record=record
  tree=reorderEdges(tree,res$tree)
  res$tree$subs=tree$edge.length
  res$tree$root.time=res$rootdate
  res$algo='TreeTime'
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
  cat(sprintf('Result of analysis using %s - model %s, clock rate %.2f, relaxation parameter %.2f, root date %.2f\n',x$algo,x$model,x$rate,x$relax,x$rootdate),...)
  invisible(x)
}

#' Plotting methods
#' @param x Output from dating
#' @param ... Additional parameters are passed on
#' @return Plot of results
#' @export
plot.resDating = function(x, ...) {
  class(x)<-'resBactDating'
  plot(x,...)
}
