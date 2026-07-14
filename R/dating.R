#' Run a dating analysis
#' @param tree Tree to date, with branch lengths in number of substitutions (cf details)
#' @param dates Dates of leaves in the tree
#' @param algo Algorithm to use, can be one of: LSD, node.dating, BactDating, treedater, TreeTime
#' @param rate If provided, force analysis to use specified value for rate, in substitutions per genome per unit of time
#' @param keepRoot Whether to keep the root as in the input
#' @param ... Passed on to dating algorithm
#'
#' @details
#' Throughout DiagnoDating, as in BactDating, the branch lengths of \code{tree} are numbers of
#' substitutions rather than substitutions per site, and the \code{rate} of a fitted model is in
#' substitutions per genome per unit of time. The \code{relax} parameter of a fitted model is on
#' the same per genome scale.
#'
#' Not all of the dating algorithms use this convention internally, and treedater in particular
#' works with rates per site. Each of them is converted, but when \code{algo='treedater'} the
#' conversion needs the length of the alignment, which should be given as \code{seqlen}, and the
#' clock model can be selected with \code{clock}. Both are passed on through \code{...}: see
#' \code{\link{runTreeDater}} for the details, which are worth reading before passing any rate
#' valued argument on to treedater.
#'
#' @return resDating object containing results of dating analysis
#' @export
#'
#' @examples
#' \dontrun{
#' #BactDating with the additive relaxed clock
#' r=runDating(tree,dates,algo='BactDating',model='arc')
#' #treedater with the same clock model, which it calls additive
#' r=runDating(tree,dates,algo='treedater',clock='additive',seqlen=10000)
#' }
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
#' @param tree Tree to date, with branch lengths in number of substitutions (cf details)
#' @param dates Dates of leaves in the tree, matched by name if named and by order otherwise
#' @param rate If provided, force analysis to use specified value for rate, in substitutions per genome per unit of time
#' @param keepRoot Whether to keep the root as in the input
#' @param seqlen Length of the alignment from which the tree was built (cf details)
#' @param clock Clock model for treedater to use, either 'strict' or 'additive'
#' @param ... Passed on to treedater::dater
#'
#' @details
#' treedater does not use the same units as the rest of DiagnoDating, so some care is needed.
#' DiagnoDating follows BactDating: the input \code{tree} has branch lengths in numbers of
#' substitutions, and \code{rate} is in substitutions per genome per unit of time. treedater
#' instead takes branch lengths in substitutions per site together with the sequence length
#' \code{s}, and reports rates in substitutions per site per unit of time.
#'
#' This function converts between the two conventions: the tree is divided by \code{seqlen}
#' on the way in, and the rate returned by treedater is multiplied by \code{seqlen} on the way
#' out. The \code{rate} argument above is therefore in DiagnoDating units (per genome), whereas
#' any rate-valued argument passed through \code{...} straight to treedater, in particular
#' \code{omega0} and \code{meanRateLimits}, is in treedater units (per site). Supplying
#' \code{seqlen} is what keeps the two consistent. If \code{seqlen} is omitted a nominal value
#' is used instead, which cancels out of the fit exactly, but which would silently rescale
#' \code{omega0} and \code{meanRateLimits} by an arbitrary factor, so passing either of those
#' without \code{seqlen} is an error.
#'
#' The clock models translate as follows. The \code{'strict'} clock of treedater is the
#' \code{'poisson'} model, which has no relaxation parameter. The \code{'additive'} clock of
#' treedater is the additive relaxed clock of Didelot et al (2021), which is exactly the
#' \code{'arc'} model of BactDating: in both, a branch of duration t under a clock of mean rate
#' mu has a negative binomial number of substitutions, with size mu*t/sigma and probability
#' 1/(1+sigma). The variance parameter that treedater calls \code{sp} is thus the same quantity
#' that BactDating calls \code{sigma} and that DiagnoDating calls \code{relax}, measured on the
#' same per genome scale, and it is carried across without rescaling. The \code{'uncorrelated'}
#' clock of treedater has no equivalent here, since its variance is not additive in the duration
#' of a branch and so cannot be written as a \code{relax} parameter.
#'
#' @return resDating object containing results of treedater analysis
#'
runTreeDater=function(tree,dates,rate=NA,keepRoot=F,seqlen=NA,clock='strict',...) {
  if (!requireNamespace('treedater',quietly=T)) stop('The treedater package is needed to use algo="treedater".')
  if (length(clock)!=1 || !clock %in% c('strict','additive')) stop(sprintf(paste0(
    'clock="%s" is not supported. Use "strict", which is the "poisson" model, or "additive", which is the "arc" model. ',
    'The "uncorrelated" clock of treedater has a variance that is not additive in the duration of a branch, ',
    'so it has no equivalent "relax" parameter.'),paste(clock,collapse=',')))

  dots=list(...)
  persite=intersect(c('omega0','meanRateLimits'),names(dots))
  if ('s' %in% names(dots)) stop('Do not pass s to treedater, use seqlen instead.')
  if (length(persite)>0 && !is.na(rate)) stop(sprintf(
    'Do not supply both rate and %s: they fix the same thing, the former per genome and the latter per site.',paste(persite,collapse=' and ')))
  if (is.na(seqlen)) {
    if (length(persite)>0) stop(sprintf(paste0(
      'seqlen is needed when passing %s on to treedater, because these are rates per site, ',
      'which cannot be interpreted without the true length of the alignment.'),paste(persite,collapse=' and ')))
    seqlen=1000 #Nominal value, which cancels out exactly as long as no per site rate is given
  }

  if (keepRoot) tre=tree else tre=unroot(tree)
  tre$edge.length=tre$edge.length/seqlen #Substitutions to substitutions per site
  sts=dates
  if (is.null(names(sts))) names(sts)=tre$tip.label #Named dates are matched by name, as in bactdate

  args=c(list(tre=tre,sts=sts,s=seqlen,clock=clock),dots)
  if (!is.na(rate)) {
    args$omega0=rate/seqlen
    args$meanRateLimits=c(1-1e-10,1+1e-10)*rate/seqlen
  }
  o=capture.output(rtd<-suppressWarnings(do.call(treedater::dater,args)))

  if (clock=='additive') {
    model='arc'
    relax=rtd$sp #The sp of treedater is the sigma of BactDating, on the same scale
    if (!is.finite(relax) || relax<=0) stop('treedater returned a relaxation parameter that is not positive, so the arc model is undefined for this fit.')
  } else {
    model='poisson'
    relax=0
  }

  res=resDating(rtd,tree,algo='treedater',model=model,rate=rtd$mean.rate*seqlen,relax=relax,rootdate=rtd$timeOfMRCA)
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
  try(suppressWarnings({
  if (keepRoot) tre=tree else tre=rtt(unroot(tree),dates)
  if (!is.na(rate)) mu=rate
  else {
    mu=ape::estimate.mu(tre,dates)
    if (mu<=0) mu=0.01
  }
  d=ape::estimate.dates(tre,dates,mu=mu)
  dt=tre
  dt$subs=dt$edge.length
  for (i in 1:nrow(dt$edge)) dt$edge.length[i]=d[dt$edge[i,2]]-d[dt$edge[i,1]]
  dt$root.time=min(d)
  res=resDating(dt,tree,algo='node.dating',model='poisson',rate=mu,relax=0,rootdate=min(d))
  return(res)}),silent = T)
  return(list(rate=NA,rootdate=NA))
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
  try(suppressWarnings({
  system(sprintf("treetime --keep-polytomies --tree /tmp/tree%d.nwk --dates /tmp/dates%d.tsv --sequence-length %d --outdir /tmp/%d %s > /dev/null",tag,tag,l,tag,opts))
  resrate=read.table(sprintf('/tmp/%d/molecular_clock.txt',tag))[1,1]*l
  restree=read.nexus(sprintf('/tmp/%d/timetree.nexus',tag))
  resrootdate=max(dates)-max(dist.nodes(restree)[Ntip(restree)+1,1:Ntip(restree)])
  restree$node.label=NULL
  trees=c(tre,restree)
  trees=.compressTipLabel(trees)
  restree=trees[[2]]
  res=resDating(restree,tree,algo='TreeTime',model='poisson',rate=resrate,relax=0,rootdate=resrootdate)
  return(res)}),silent = T)
  return(list(rate=NA,rootdate=NA))
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
