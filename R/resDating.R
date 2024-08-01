#' Constructor for the class resDating
#'
#' @param dt Dated phylogeny
#' @param phy Undated phylogeny (may be rooted or not)
#' @param algo Algorithm used for dating
#' @param model Model used for dating
#' @param rate Clock rate used for dating
#' @param relax Relaxation parameter used for dating
#' @param rootdate Date of the root inferred by dating
#'
#' @return Object of class resDating
#' @export
#'
resDating = function(dt, phy, algo='Unknown', model='poisson', rate=10, relax=0, rootdate=NA)
{
  phy=unroot(phy)
  r=list()
  r$rate=rate
  r$relax=relax
  r$tree=dt
  if (is.na(rootdate)) r$rootdate=dt$root.time else r$rootdate=rootdate
  r$model=model
  r$algo=algo
  class(r) <- 'resDating'
  subs=rep(NA,length(dt$edge.length))
  matphy=matrixChildren(phy)
  matdt=matrixChildren(dt)
  for (i in 1:nrow(dt$edge)) {
    prop=1
    if (dt$edge[i,1]==Ntip(dt)+1) #Branch connected to root
      prop=dt$edge.length[i]/sum(dt$edge.length[which(dt$edge[,1]==Ntip(dt)+1)])
    cdt=dt$tip.label[matdt[dt$edge[i,2],]]
    found=NA
    for (j in 1:nrow(phy$edge)) {
      cphy=phy$tip.label[matphy[phy$edge[j,2],]]
      if (setequal(cphy,cdt) || setequal(cphy,setdiff(dt$tip.label,cdt))) {
        subs[i]=phy$edge.length[j]*prop
      }
    }
  }
  r$tree$subs=subs
  if (abs(sum(r$tree$subs)-sum(phy$edge.length))>0.01) warning('Incorrect number of subs.')
  if (any(is.na(r$tree$subs))) warning('NAs in subs, maybe change in topology.')
  return(r)
}

#' Print function for resDating objects
#' @param x output from dating
#' @param ... Passed on to cat
#' @return Print out details of dating results
#' @export
print.resDating <- function(x, ...)
{
  stopifnot(inherits(x, "resDating"))
  cat(sprintf('Result from %s, model %s, clock rate %.2f, relaxation parameter %.2f, root date %.2f\n',x$algo,x$model,x$rate,x$relax,x$rootdate),...)
  invisible(x)
}

#' Plotting methods
#' @param x Output from dating
#' @param ... Additional parameters are passed on
#' @return Plot of results
#' @export
plot.resDating = function(x, ...) {
  stopifnot(inherits(x, "resDating"))
  plot(x$tree,...)
  axisPhylo(backward = F)
}
