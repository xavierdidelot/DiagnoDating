#' Test to determine whether there are lineages with unbalanced dates
#' For each node, perform a Kruskal-Wallis test for the dates on the left and on the right
#'
#' @param tree Tree to date
#' @param dates Dates of leaves in the tree
#'
#' @return p-values
#' @export
#'
testConfounding = function(tree,dates) {
  children=matrixChildren(tree)
  pvals=rep(NA,Nnode(tree))
  for (node in 1:Nnode(tree)) {
    w=tree$edge[which(tree$edge[,1]==node+Ntip(tree)),2]
    dates1=dates[children[w[1],]]
    dates1=dates1[!is.na(dates1)]
    dates2=dates[children[w[2],]]
    dates2=dates2[!is.na(dates2)]
    res=kruskal.test(c(dates1,dates2),c(rep(1,length(dates1)),rep(2,length(dates2))))
    pvals[node]=res$p.value
    if (res$p.value<1e-5) {
      print(dates1)
      print(dates2)
    }
  }

  return(pvals)
}
