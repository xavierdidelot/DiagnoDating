#Test probabilities
context("Test plots")

test_that("All plots work.", {
  set.seed(0)
  dates=runif(10,2000,2020)
  dt=simcoaltree(dates)
  phy=simobsphy(dt,model='arc',mu=10,sigma=5)
  expect_silent(r<-resDating(dt,phy,model='arc',rate=10,relax=5))
  expect_is(r,'resDating')
  expect_silent(plot(r))
  expect_silent(plotResid(r))
  expect_silent(plotResid(r,sub=1,xlab='',ylab='',main=''))
  expect_silent(plotResid(r,sub=2,xlab='',ylab='',main=''))
  expect_silent(plotResid(r,sub=3,xlab='',ylab='',main=''))
  expect_silent(plotResid(r,sub=4,xlab='',ylab='',main=''))
  expect_silent(plotProbBranches(r))
  expect_silent(plotProbBranches(r,color=F))
  expect_silent(plotProbBranches(r,sub=1,xlab='',ylab='',main=''))
  expect_silent(plotProbBranches(r,sub=2,main=''))
})

