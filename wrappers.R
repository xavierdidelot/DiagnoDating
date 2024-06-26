runTreeDater=function(tree,dates) {
  tre=tree
  l=1e6
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

runLSD=function(tree,dates) {
  tre=tree
  l=1e6
  tre$edge.length=tre$edge.length/l
  tre$tip.label=1:Ntip(tre)
  sts=dates
  names(sts)=1:Ntip(tre)
  write.tree(tre,'/tmp/tree.nwk')
  write.table(sts,'/tmp/dates.csv',quote = F,col.names=length(sts))
  system('R << ')
  o=capture.output(result <- Rlsd2::lsd2(inputTree='/tmp/tree.nwk', inputDate='/tmp/dates.csv',outFile = '/tmp/result', seqLen=l))
  res=list()
  res$rate=result$rate*l
  res$tmrca=result$tMRCA
  return(res)
}