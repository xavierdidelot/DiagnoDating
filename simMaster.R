#' Simulation with structure using Master
#' @param demes Number of demes
#' @param coalrate Within-deme coalescent rate
#' @param migrate Migration rate
#' @param samplingDates A list containing for each population a vector of sampling dates
#' @return A simulated dated phylogeny
#' @export
simMaster = function(demes=3,coalrate=1,migrate=0.1,NeFunImp,samplingDates)
{
  #Generate XML file for Master
  f=file('/tmp/master2.xml',open='w')
  writeLines("<beast version='2.0' namespace='master:master.model:master.steppers:master.conditions:master.postprocessors:master.outputs'>",f)
  writeLines("<run spec='InheritanceTrajectory' verbosity='2'>",f)
  writeLines("<model spec='Model'>",f)
  writeLines(sprintf("<populationType spec='PopulationType' typeName='L' id='L' dim='%d'/>",demes),f)

  writeLines("<reactionGroup spec='ReactionGroup' reactionGroupName='Coalescence'>",f)
  for (i in 1:demes)
      writeLines(sprintf("<reaction spec='Reaction' rate='%f'>2L[%d]:1 -> L[%d]:1</reaction>",coalrate,i-1,i-1),f)
  writeLines("</reactionGroup>",f)
  writeLines("<reactionGroup spec='ReactionGroup' reactionGroupName='Migration'>",f)
  for (i in 1:demes) for (j in 1:demes) if (i!=j)
    writeLines(sprintf("<reaction spec='Reaction' rate='%f'>L[%d] -> L[%d]</reaction>",migrate,i-1,j-1),f)
  writeLines("</reactionGroup>",f)
  writeLines("</model>",f)
  writeLines("<initialState spec='InitState'>",f)

  dates=unlist(samplingDates)
  locs=c()
  for (i in 1:length(samplingDates)) locs=c(locs,rep(i-1,length(samplingDates[[i]])))
  for (i in 1:length(dates)) {
    l=sprintf('\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f"><population spec="Population" type="@L" location="%d"/></lineageSeedMultiple>\n',dates[i],locs[i])
    writeLines(l,f)
  }
  writeLines('</initialState>',f)
  writeLines('  <lineageEndCondition spec="LineageEndCondition" nLineages="1"/>',f)
  writeLines('  <output spec="NexusOutput" fileName="/tmp/simu_master.tree" collapseSingleChildNodes="false" reverseTime="true"/>',f)
               writeLines('  </run>',f)
  writeLines('</beast>',f)
  close(f)

  #Generate tree with Master
  system(sprintf('beast2 -seed %d /tmp/master2.xml > /dev/null 2> /dev/null',round(runif(1,0,1e6))))
  t=treeio::read.beast('/tmp/simu_master.tree')
  phy=as.phylo(t)
  phy=collapse.singles(phy)#remove migration nodes
  return(phy)
}
