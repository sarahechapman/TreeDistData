context("Loading data")

test_that("Data dimensions are correct", {
  message('D-4')
  data('tdMethods', package = 'TreeDistData')
  nMetrics <- length(tdMethods)

  message('D-8')
  LengthWithout <- function(...) nMetrics - length(list(...))

  message('D-7td')
  data("sevenTipDistances", package='TreeDistData')
  expect_equal(LengthWithout('mafi'), length(sevenTipDistances))
  expect_equal(c(945L, 945L), dim(sevenTipDistances[[1]]))

  message('D-dd25')
  data("distanceDistribution25", package='TreeDistData')
  expect_equal(c(LengthWithout('mafi'), 10000L),
               dim(distanceDistribution25))
  expect_equal(dim(distanceDistribution50), dim(distanceDistribution25))

  noMafiNnit <- tdMethods[!tdMethods %in% c('mafi', 'nni_t')]

  message('D-rtd')
  data("randomTreeDistances", package='TreeDistData')
  nLeafMeasurements <- 197L
  expect_true(all(noMafiNnit %in% dimnames(randomTreeDistances)[[1]]))
  expect_equal(c(LengthWithout('mafi', 'nni_t'), 13L, nLeafMeasurements),
               dim(randomTreeDistances))

  message('D-bull')
  expect_true(all(noMafiNnit %in% names(bullseyeDistances[[1]])))
  expect_true(all(noMafiNnit %in% dimnames(bullseyeMorphScores[[1]])[[2]]))
  #TODO update other datasets for all 20 metrics?
  message('D-one')
})
