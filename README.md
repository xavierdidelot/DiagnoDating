
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

Work on the confounding effect of genetic structure on molecular dating.

# Done so far

- Reproduce confounding effect described by Murray2016 on similar
  coalescent simulation and using BactDating and treedater, cf
  `confounding.Rmd`
- Ruled out possible problem in phylogenetics since working directly
  with trees, not genetic data
- Ruled out possible problem with root identification (ie consider root
  known)
- Ruled out possible problem with multifurcation at the root (relaxed
  this in simulation)
- More general version in `confounding-simstructure.Rmd` in which
  sampling dates do not have to be identical, and population sizes do
  not have to be constant. Simulation code is in `simStructure.R`, uses
  some of mlesky code for simulation of non-homogeneous Poisson
  coalescent process.
- Difficult to get confusion using DetectImports simulations, cf
  `confounding-detectimports.Rmd`. May be because confusion only happens
  when sampling dates are biased, whereas DetectImports assumes sampling
  according to relative prevalence of populations
- Mendel test will flag an issue even if there is no structure and no
  problem with dating, cf `struture-test.Rmd`
- When there is no confounding issue but the temporal signal is weak, eg
  because all sampling dates are the same or very similar, see
  `run-nosignal.Rmd`. Estimate incorrectly high clock rate. But the
  permutation test shows this is not to be trusted.

# Still to do

- Propose new test to detect when the confounding problem happens
- Propose a solution to when the confounding problem happens, maybe via
  removal of some leaves until test satisfied
- Is permutation test useful?
- Is clustered permutation test useful?
- Comparison with run in which all leaves have same date. Attractive as
  only need to run once. Use DIC or something else?
- Use better Bayesian model comparison methods?
- Use Bayesian model criticism approach? Maybe using posterior
  predictive p-values?
- Effect of relaxed clock models, especially ARC.
- Effect of priors, especially on tree (constant effective population
  size as in Murray2016) but also others.
- Apply to real data. May be able to reuse data from Holden2013, cf Fig
  S4 in Murray2016.

# Some ideas

- Fit with all dates equal using modified BactDating. Need to add move
  to scale up tree and down rate simultaneously, and vice versa.
  Consider prior effect. Given this fit, do posterior predictive test to
  see if there is a temporal signal.
- Joint model with and without temporal signal, to estimate Bayes Factor
  directly with rjMCMC. Current moves on sampling times could be
  starting point.
- Show confounding only happens when we have structure plus uneven
  sampling between structure component. Develop test for this based on
  undated phylgeny plus sampling dates. Bit like treeBreaker with
  continuous phenotype representing sampling dates.
- Can a relaxed clock with high relaxation parameter equate a lack of
  temporal signal? What happens in ARC if omega is very high?
