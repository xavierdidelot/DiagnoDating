#Make matrix of children
matrixChildren = function(tree) {
  parent=NA
  parent[tree$edge[,2]]=tree$edge[,1]
  children=matrix(F,Ntip(tree)+Nnode(tree),Ntip(tree))
  for (i in 1:Ntip(tree)) {
    children[i,i]=T
    cur=i
    while (cur!=Ntip(tree)+1) {
      cur=parent[cur]
      children[cur,i]=T
    }
  }
  return(children)
}


#Reorder the edges of t1 so that they are the same as in t2
reorderEdges=function(t1,t2) {
  mat1=matrixChildren(t1)
  mat2=matrixChildren(t2)
  inds=rep(NA,nrow(t2$edge))
  for (i in 1:nrow(t2$edge)) {
    c2=t2$tip.label[mat2[t2$edge[i,2],]]
    found=NA
    for (j in 1:nrow(t1$edge)) {
      c1=t1$tip.label[mat1[t1$edge[j,2],]]
      if (setequal(c1,c2)) {
        found=j
        break
      }
    }
    inds[i]=found
  }
  t1$edge=t1$edge[inds,]
  t1$edge.length=t1$edge.length[inds]
  return(t1)
}