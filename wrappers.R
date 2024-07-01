runTreeDater=function(tree,dates) {
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l)))
  res=list()
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
  class(res)<-'resBactDating'
  return(res)
}

runLSD=function(tree,dates,tag=0) {
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  write.tree(tre,sprintf('/tmp/tree%d.nwk',tag))
  write.table(sts,sprintf('/tmp/dates%d.csv',tag),quote = F,col.names=length(sts))
  o=capture.output(result <- Rlsd2::lsd2(inputTree=sprintf('/tmp/tree%d.nwk',tag), inputDate=sprintf('/tmp/dates%d.csv',tag),outFile = sprintf('/tmp/result%d',tag), seqLen=l))
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
  res$tree=multi2di(res$tree)
  class(res)<-'resBactDating'
  return(res)
}