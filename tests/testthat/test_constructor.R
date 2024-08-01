#Test the constructor function
context("Test the constructor function")

test_that("Constructor works.", {
  set.seed(0)
  dt=simcoaltree(runif(10,2000,2020))
  p=simobsphy(dt)
  p=root(p,"1")
  expect_silent(res<-resDating(dt,p))
  expect_is(res,'resDating')
})

