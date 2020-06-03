context("Rootedness")

test_that("Root position doesn't influence tree score", {

  library('TreeTools')
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)

  trees <- lapply(list(bal8, RootTree(bal8, 't4'),
                       pec8, RootTree(pec8, 't4')), UnrootTree)

  expect_equal(AllDists(trees[[1]], trees[[2]]),
               AllDists(trees[[1]], trees[[1]]),
               tolerance = 1e-7)

  expect_equal(AllDists(trees[[3]], trees[[4]]),
               AllDists(trees[[3]], trees[[3]]),
               tolerance = 1e-7)

  expect_equal(AllDists(trees[[1]], trees[[3]]),
               AllDists(trees[[2]], trees[[4]]),
               tolerance = 1e-7)

  expect_equal(AllDists(trees[[1]], trees[[3]]),
               AllDists(trees[[1]], trees[[4]]),
               tolerance = 1e-7)

  expect_equal(AllDists(trees[[1]], trees[[4]]),
               AllDists(trees[[2]], trees[[4]]),
               tolerance = 1e-7)

})
