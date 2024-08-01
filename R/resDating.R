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
resDating = function(dt, phy, algo='Unknown', model='poisson', rate=10, relax=NA, rootdate=NA)
{
  r=list()
  r$rate=rate
  r$relax=relax
  r$tree=dt
  if (is.na(rootdate)) r$rootdate=dt$root.time else r$rootdate=rootdate
  r$model=model
  r$algo='algo'
  #r$tree$subs=tree$edge.length#TODO
  class(r) <- 'resDating'
  return(r)
}

#' Print function for resDating objects
#' @param x output from bactdate
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
  class(x)<-'resBactDating'
  plot(x,...)
}
