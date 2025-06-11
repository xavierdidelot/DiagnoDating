#Test probabilities
context("Test probabilities")

test_that("Probabilities agree with simulation.", {
  set.seed(0)
  samplingDates=list(2023:2025,2023:2025,2023:2025)
  t=simStructure(samplingDates = samplingDates)

  #Poisson
  p=simobsphy(t,model='poisson',mu=10)
  r=list(model='poisson',rate=10,relax=0,tree=t)
  r$tree$subs=p$edge.length
  class(r) <- 'resDating'
  expect_equal(p$prob,sum(calcLikBranches(r,log=T)))

  #ARC
  p=simobsphy(t,model='arc',mu=10,sigma=2)
  r=list(model='arc',rate=10,relax=2,tree=t)
  r$tree$subs=p$edge.length
  class(r) <- 'resDating'
  expect_equal(p$prob,sum(calcLikBranches(r,log=T)))

  #Strict Gamma
  p=simobsphy(t,model='strictgamma',mu=10)
  r=list(model='strictgamma',rate=10,relax=0,tree=t)
  r$tree$subs=p$edge.length
  class(r) <- 'resDating'
  expect_equal(p$prob,sum(calcLikBranches(r,log=T)))

  #cARC
  p=simobsphy(t,model='carc',mu=10,sigma=2)
  r=list(model='carc',rate=10,relax=2,tree=t)
  r$tree$subs=p$edge.length
  class(r) <- 'resDating'
  expect_equal(p$prob,sum(calcLikBranches(r,log=T)))

})

