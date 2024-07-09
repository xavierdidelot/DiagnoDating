#Test phylo utils
context("Test phylo utils")

test_that("Reordering tree edges according to itself does not make a change.", {
  t=rtree(10)
  expect_equal(t$edge,reorderEdges(t,t)$edge)
})

test_that("Root has all leaves as children.", {
  t=rtree(10)
  expect_silent(m<-matrixChildren(t))
  expect_equal(sum(m[Ntip(t)+1,]),Ntip(t))
})

