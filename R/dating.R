#' Run a dating analysis
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param algo Algorithm to use, can be one of: LSD, node.dating, BactDating, treedater, TreeTime
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on to dating algorithm
#'
#' @return resDating object containing results of dating analysis
#' @export
#'
runDating=function(tree,dates,algo='BactDating',rate=NA,keepRoot=F,...) {
  r=NULL
  if (!is.rooted(tree) && keepRoot) stop('Need rooted tree as input for keepRoot option')
  if (algo=='1' || algo=='LSD') r=runLSD(tree,dates,rate,keepRoot,...)
  if (algo=='2' || algo=='node.dating') r=runNodeDating(tree,dates,rate,keepRoot,...)
  if (algo=='3' || algo=='BactDating') r=runBactDating(tree,dates,rate,keepRoot,...)
  if (algo=='4' || algo=='treedater') r=runTreeDater(tree,dates,rate,keepRoot,...)
  if (algo=='5' || algo=='TreeTime') r=runTreeTime(tree,dates,rate,keepRoot,...)
  if (is.null(r)) stop('Unknown algorithm')
  return(r)
}

#' Date a tree using BactDating
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on to BactDating::bactdate
#'
#' @return resDating object containing results of BactDating analysis
#'
runBactDating=function(tree,dates,rate=NA,keepRoot=F,...) {
  def_args=list(tree=tree,date=dates,model='poisson')
  if (!is.na(rate)) def_args=c(def_args,initMu=rate,updateMu=F)
  if (keepRoot) def_args=c(def_args,updateRoot=F)
  cl=as.list(match.call())[-(1:5)]
  args=c(cl,def_args[!names(def_args) %in% names(cl)])
  r=do.call("bactdate",args)
  r$algo='BactDating'
  r$model=args$model
  attributes(r$tree)$order<-NULL
  attributes(r$inputtree)$order<-NULL
  v=r$record[,'mu']
  v=v[(1+length(v)/2):length(v)]
  r$rate=mean(v)
  v=r$record[,'sigma']
  v=v[(1+length(v)/2):length(v)]
  r$relax=mean(v)
  class(r)<-'resDating'
  r$resid=calcResiduals(r)
  return(r)
}

#' Extract from BactDating output (r) the tree corresponding to a specific iteration in (w)
#' @param r Output from BactDating
#' @param w Iteration to extract (default is the last iteration)
#'
#' @return BactDating object corresponding to specified iteration
#' @export
#'
takeSample=function(r,w=nrow(r$record)) {
  tree = r$tree
  rowRec = r$record[w,]
  s=Ntip(tree)+Nnode(tree)
  i2=1:nrow(tree$edge)
  tree$edge[i2,1]=r$record[w,s+tree$edge[i2,2]]
  tree$edge.length[i2] = rowRec[tree$edge[i2, 2]] - rowRec[tree$edge[i2, 1]]
  tree$subs[i2] = rowRec[s*2+tree$edge[i2,2]]
  rmod=r
  rmod$tree=tree
  rmod$rootdate=r$record[w,Ntip(tree)+1]
  rmod$rate=r$record[w,'mu']
  rmod$relax=r$record[w,'sigma']
  rmod$resid=calcResiduals(rmod)
  return(rmod)
}

#' Date a tree using treedater
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on to treedater::dater
#'
#' @return resDating object containing results of treedater analysis
#'
runTreeDater=function(tree,dates,rate=NA,keepRoot=F,...) {
  if (keepRoot) tre=tree else tre=unroot(tree)
  l=1000#round(sum(tre$edge.length)*100)
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tre$tip.label
  if (is.na(rate)) o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l,...)))
  else o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l,omega0=rate/l,meanRateLimits=c(1-1e-10,1+1e-10)*rate/l,...)))
  res=resDating(rtd,tree,algo='treedater',model='poisson',rate=rtd$mean.rate*l,relax=0,rootdate=rtd$timeOfMRCA)
  return(res)
}

#' Date a tree using LSD
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Ignored for the time being
#'
#' @return resDating object containing results of LSD analysis
#'
runLSD=function(tree,dates,rate=NA,keepRoot=F,...) {
  tag=round(runif(1,1,1e8))
  if (keepRoot) tre=tree else tre=unroot(tree)
  l=round(sum(tre$edge.length)*100)
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tre$tip.label
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.csv',tag),quote = F,col.names=length(sts))
  if (keepRoot) opts='' else opts='-r a'
  if (is.na(rate)) system(sprintf("lsd2 -i /tmp/tree%d.nwk -d /tmp/dates%d.csv -s %d -l -1 %s > /dev/null",tag,tag,l,opts))
  else system(sprintf("echo %f > /tmp/rate%d;lsd2 -i /tmp/tree%d.nwk -d /tmp/dates%d.csv -s %d -l -1 %s -w /tmp/rate%d > /dev/null",rate/l,tag,tag,tag,l,opts,tag))
  lines=readLines(sprintf('/tmp/tree%d.nwk.result',tag))
  lines=lines[grep('tMRCA',lines)]
  lines=as.numeric(unlist(strsplit(lines, "[ ,]"))[c(3,6)])
  rtd=read.nexus(sprintf('/tmp/tree%d.nwk.result.date.nexus',tag))
  trees=c(tre,rtd)
  trees=.compressTipLabel(trees)
  rtd=trees[[2]]
  res=resDating(rtd,tree,algo='LSD',model='poisson',rate=lines[1]*l,relax=0,rootdate=lines[2])
  return(res)
}

#' Date a tree using node.dating
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on to ape::node.dating
#'
#' @return resDating object containing results of node.dating analysis
#'
runNodeDating=function(tree,dates,rate=NA,keepRoot=F,...) {
  if (keepRoot) tre=tree else try(suppressWarnings(tre<-rtt(unroot(tree),dates)),silent=T)
  if (!is.na(rate)) mu=rate
  else {
    try(suppressWarnings(mu<-ape::estimate.mu(tre,dates)),silent=T)
    if (mu<=0) mu=0.01
  }
  suppressWarnings(d<-ape::estimate.dates(tre,dates,mu=mu))
  dt=tre
  dt$subs=dt$edge.length
  for (i in 1:nrow(dt$edge)) dt$edge.length[i]=d[dt$edge[i,2]]-d[dt$edge[i,1]]
  dt$root.time=min(d)
  res=resDating(dt,tree,algo='node.dating',model='poisson',rate=mu,relax=0,rootdate=min(d))
  return(res)
}

#' Date a tree using TreeTime
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#' @param rate If provided, force analysis to use specified value for rate
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on
#'
#' @return resDating object containing results of TreeTime analysis
#'
runTreeTime=function(tree,dates,rate=NA,keepRoot=F,...) {
  tag=round(runif(1,1,1e8))
  if (keepRoot) tre=tree else tre=unroot(tree)
  l=round(sum(tre$edge.length)*100)
  tre$edge.length=tre$edge.length/l
  sts=dates
  names(sts)=tre$tip.label
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.tsv',tag),quote = F,col.names='strain\tdate',sep='\t')
  if (keepRoot) opts='--keep-root' else opts=''
  if (!is.na(rate)) opts=sprintf('%s --clock-rate %f',opts,rate/l)
  system(sprintf("treetime --keep-polytomies --tree /tmp/tree%d.nwk --dates /tmp/dates%d.tsv --sequence-length %d --outdir /tmp/%d %s > /dev/null",tag,tag,l,tag,opts))
  resrate=read.table(sprintf('/tmp/%d/molecular_clock.txt',tag))[1,1]*l
  restree=read.nexus(sprintf('/tmp/%d/timetree.nexus',tag))
  resrootdate=max(dates)-max(dist.nodes(restree)[Ntip(restree)+1,1:Ntip(restree)])
  restree$node.label=NULL
  trees=c(tre,restree)
  trees=.compressTipLabel(trees)
  restree=trees[[2]]
  res=resDating(restree,tree,algo='TreeTime',model='poisson',rate=resrate,relax=0,rootdate=resrootdate)
  return(res)
}

#' Estimation of the coalescent rate alpha
#'
#' @param tree Dated tree
#' @param sampling Whether to sample from the distribution or return the mean
#'
#' @return Estimated value of the coalescent rate alpha
#' @export
#'
estimAlpha=function(tree,sampling=F)
{
  n=Ntip(tree)
  nr=2*n-1
  d=dist.nodes(tree)[n+1,]
  s=sort(d,decreasing = T, index.return = TRUE)
  k=cumsum(2*(s$ix<=n)-1)
  difs=s$x[1:(nr-1)]-s$x[2:nr]
  su=sum(k[1:(nr-1)]*(k[1:(nr-1)]-1)*difs)
  #Using inverse-gamma prior
  if (sampling==F) alpha=1/((n-1)*2/su)
  else alpha=1/rgamma(1,shape=n-1,scale=2/su)
  return(alpha)
}
