# New test to determine whether there are lineages with unbalanced dates
testConfounding = function(tree,dates) {
  #Make vector of parents
  parent=NA
  parent[tree$edge[,2]]=tree$edge[,1]
  
  #Make matrix of children
  children=matrix(NA,Ntip(tree)+Nnode(tree),Ntip(tree))
  for (i in 1:Ntip(tree)) {
    children[i,1]=i
    cur=i
    while (cur!=Ntip(tree)+1) {
      cur=parent[cur]
      w=match(NA,children[cur,])
      children[cur,w]=i
    }
  }
  
  #For each node, perform a Kruskal-Wallis test for the dates on the left and on the right
  pvals=rep(NA,Nnode(tree))
  for (node in 1:Nnode(tree)) {
    w=tree$edge[which(tree$edge[,1]==node+Ntip(tree)),2]
    dates1=dates[children[w[1],]]
    dates1=dates1[!is.na(dates1)]
    dates2=dates[children[w[2],]]
    dates2=dates2[!is.na(dates2)]
    res=kruskal.test(c(dates1,dates2),c(rep(1,length(dates1)),rep(2,length(dates2))))
    pvals[node]=res$p.value
  }
  
  return(pvals)
}