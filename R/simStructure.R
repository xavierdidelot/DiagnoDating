#' Simulation with structure
#' @param popStarts Vector of dates when populations start
#' @param globalNeg Value of Ne*g for the global population
#' @param NeFunPop Function determining the demographic trajectory of populations, default is CaveDive model with N=2 and h=0.2
#' @param samplingDates A list containing for each population a vector of sampling dates
#' @return A simulated dated phylogeny
#' @export
simStructure = function(popStarts=c(2020,2021,2022),globalNeg=1,NeFunPop,samplingDates)
{
  #Create Neg functions for each population
  npop=length(popStarts)
  if (missing(NeFunPop)) NeFunPop=function(t,texp) {(t>texp)*2*(t-texp)^2/(0.2^2+(t-texp)^2)}
  NeFun=list()
  for (i in 1:npop)
    NeFun[[i]]=function(t) {NeFunPop(t,popStarts[i])}
  
  #Simulate population trees
  popTrees=list(NA,npop)
  toadd=rep(NA,npop)
  for (i in 1:npop) {
    while (is.na(toadd[i])||toadd[i]<0) {
      popTrees[[i]]=mlesky::simCoal(samplingDates[[i]],NeFun[[i]],1e-2)
      toadd[i]=popTrees[[i]]$root.time-popStarts[i]}
  }
  
  #Case without structure
  if (npop==1) {
    t=popTrees[[1]]
    return(t)
  }
  
  #Simulate global tree
  t=mlesky::simCoal(popStarts,function(t){return(globalNeg)})
  t$tip.label=sprintf('G%d',1:Ntip(t))
  
  #Paste trees together
  a=0
  subpops=rep(NA,npop-1)
  for (i in 1:npop) {
    if (i>1) subpops[i-1]=a+which(samplingDates[[i]]==min(samplingDates[[i]]))[1]
    t2=popTrees[[i]]
    t2$tip.label=as.numeric(a+(1:Ntip(t2)))
    w=which(t$tip.label==sprintf('G%d',i))
    w2=which(t$edge[,2]==w)
    t$edge.length[w2]=t$edge.length[w2]+toadd[i]
    t=bind.tree(t,t2,where=w,position=0)
    a=a+Ntip(t2)
  }
  return(t)
}
