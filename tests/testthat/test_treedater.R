#Test the treedater interface, and especially the conversion of units
context("Test treedater interface")

test_that("treedater runs under both clock models and reports the right parameters", {
  skip_if_not_installed('treedater')
  set.seed(0)
  L=10000
  phy=simcoaltree(dates=2024+runif(30,0,2),alpha=2)
  obs=simobsphy(phy,mu=10,sigma=0.5,model='arc')#Branch lengths in substitutions
  dates=leafDates(phy)

  #Strict clock is the Poisson model, and has no relaxation parameter
  expect_silent(r<-runDating(obs,dates,algo='treedater',seqlen=L))
  expect_is(r,'resDating')
  expect_equal(r$model,'poisson')
  expect_equal(r$relax,0)
  expect_true(all(is.finite(r$resid)))

  #Additive clock is the arc model, and its sp is the relaxation parameter
  expect_silent(r<-runDating(obs,dates,algo='treedater',seqlen=L,clock='additive'))
  expect_equal(r$model,'arc')
  expect_true(r$relax>0)
  expect_true(all(is.finite(r$resid)))#Would be NaN if relax was left at zero
  expect_is(testResid(r),'htest')
})

test_that("treedater results do not depend on the nominal sequence length", {
  skip_if_not_installed('treedater')
  set.seed(0)
  phy=simcoaltree(dates=2024+runif(30,0,2),alpha=2)
  obs=simobsphy(phy,mu=10,sigma=0.5,model='arc')
  dates=leafDates(phy)

  #The rate is per genome and the relaxation parameter is per genome, so neither of them should
  #move when the tree is handed to treedater cut into a different number of pieces
  r1=runDating(obs,dates,algo='treedater',seqlen=1000,clock='additive')
  r2=runDating(obs,dates,algo='treedater',seqlen=50000,clock='additive')
  expect_equal(r1$rate,r2$rate)
  expect_equal(r1$relax,r2$relax)

  #A rate given to runDating is per genome, and should come back unchanged
  r3=runDating(obs,dates,algo='treedater',seqlen=10000,rate=10)
  expect_equal(r3$rate,10,tolerance=1e-4)
})

test_that("treedater dates are matched by name when they are named", {
  skip_if_not_installed('treedater')
  set.seed(0)
  phy=simcoaltree(dates=2024+runif(20,0,2),alpha=2)
  obs=simobsphy(phy,mu=10,sigma=0.5,model='arc')
  dates=leafDates(phy)
  names(dates)=obs$tip.label

  r1=runDating(obs,dates,algo='treedater',seqlen=10000)
  r2=runDating(obs,dates[sample(length(dates))],algo='treedater',seqlen=10000)
  expect_equal(r1$rootdate,r2$rootdate)
})

test_that("treedater arguments that are per site are refused without a sequence length", {
  skip_if_not_installed('treedater')
  set.seed(0)
  phy=simcoaltree(dates=2024+runif(20,0,2),alpha=2)
  obs=simobsphy(phy,mu=10,sigma=0.5,model='arc')
  dates=leafDates(phy)

  expect_error(runDating(obs,dates,algo='treedater',meanRateLimits=c(1e-3,2e-3)),'seqlen')
  expect_error(runDating(obs,dates,algo='treedater',omega0=1e-3),'seqlen')
  #The uncorrelated clock of treedater has no equivalent relaxation parameter
  expect_error(runDating(obs,dates,algo='treedater',seqlen=1e4,clock='uncorrelated'),'not supported')
  #The rate is per genome and omega0 is per site, so giving both is contradictory
  expect_error(runDating(obs,dates,algo='treedater',seqlen=1e4,rate=10,omega0=1e-3),'both')
})
