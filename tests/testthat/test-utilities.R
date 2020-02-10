context('utilities.R')

test_that('Pairwise distances calculated correctly', {
  nTrees <- 6L
  nTip <- 16L

  message('U-7')
  set.seed(0)
  trees <- lapply(rep(nTip, nTrees), ape::rtree, br=NULL)
  trees[[1]] <- TreeTools::BalancedTree(nTip)
  trees[[nTrees - 1L]] <- TreeTools::PectinateTree(nTip)
  class(trees) <- 'multiPhylo'

  message('U-14')
  dists <- PairwiseDistances(trees, phangorn::RF.dist)
  expect_equivalent(phangorn::RF.dist(trees), dists)

  message('U-18')
  lapply(CompareAllTrees(trees), function (dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })


  message('U-24')
  sprWalk <- vector('list', nTrees)
  sprWalk[[1]] <- lastTree <- TreeTools::PectinateTree(nTip)

  message('U-28')
  for (i in seq.int(2, nTrees)) {
    sprWalk[[i]] <- lastTree <- TreeSearch::SPR(lastTree)
  }

  message('U-33')
  lapply(CompareAllTrees(sprWalk), function(dist) {
    expect_equal(c(nTrees, nTrees), dim(as.matrix(dist)))
  })
  message('U-##')
})
