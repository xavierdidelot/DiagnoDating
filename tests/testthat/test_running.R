#Test all is running without error
context("Test running without error")

test_that("Simulate with simStructure and inference with all methods", {
  set.seed(0)
  samplingDates=list(2023:2025,2023:2025,2023:2025)
  t=simStructure(samplingDates = samplingDates)
  p=simobsphy(t)
  expect_silent(res<-runDating(p,unlist(samplingDates),algo='BactDating',nbIts=10))
  expect_is(res,'resDating')
  expect_is(capture_output(print(res)),'character')
#  expect_silent(res<-runDating(p,unlist(samplingDates),algo='LSD'))
#  expect_is(res,'resDating')
#  expect_is(capture_output(print(res)),'character')
  expect_silent(res<-runDating(p,unlist(samplingDates),algo='node.dating'))
  expect_is(res,'resDating')
  expect_is(capture_output(print(res)),'character')
  expect_silent(res<-runDating(p,unlist(samplingDates),algo='treedater'))
  expect_is(res,'resDating')
  expect_is(capture_output(print(res)),'character')
  expect_silent(res<-runDating(p,unlist(samplingDates),algo='TreeTime'))
  expect_is(res,'resDating')
  expect_is(capture_output(print(res)),'character')
  expect_silent(plotProbBranches(res))
  expect_silent(plotResid(res))
  expect_silent(t<-testResid(res))
  expect_is(t,'htest')
})

test_that("Simulate with simMaster.", {
  set.seed(0)
  samplingDates=list(2023:2025,2023:2025,2023:2025)
  t=simMaster(samplingDates = samplingDates)
  p=simobsphy(t)
  expect_silent(res<-runDating(p,unlist(samplingDates),nbIts=10))
  expect_is(res,'resDating')
  expect_is(capture_output(print(res)),'character')
})

test_that("testImbalance works.", {
  set.seed(0)
  p=rtree(10)
  d=2001:2010
  expect_silent(res<-testImbalance(p,d))
  expect_is(res,'numeric')
  expect_equal(length(res),length(d)-1)
})

