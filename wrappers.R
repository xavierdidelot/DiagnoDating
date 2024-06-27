runTreeDater=function(tree,dates) {
  tre=tree
  l=max(tree$edge.length)*1000
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  o=capture.output(rtd<-suppressWarnings(treedater::dater(tre,sts,s=l)))
  res=list()
  res$rate=rtd$mean.rate*l
  res$tmrca=rtd$timeOfMRCA
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
  res$rate=result$rate*l
  res$tmrca=result$tMRCA
  return(res)
}