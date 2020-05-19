library('TreeDist')
library('TreeTools')
library('TreeSearch')
library('TreeDistData')
data('bullseyeTrees', package = 'TreeDistData')

Successes <- function (Distance, trees) {
  nTip <- NTip(trees[[1]])
  sum(vapply(trees, function (tr) {
    GenerateSet <- function (Rearrange, n) {
      structure(c(lapply(rep(0, n), Rearrange, tr = tr), list(tr)),
                class = 'multiPhylo')
    }
    llis <- GenerateSet(function (O, tr)
      LeafLabelInterchange(tr, ceiling(nTip * .4)), 7L)
    nnis <- GenerateSet(function (O, tr) NNI(NNI(NNI(NNI(tr)))), 7L)
    sprs <- GenerateSet(function (O, tr) SPR(SPR(SPR(SPR(tr)))), 7L)
    tbrs <- GenerateSet(function (O, tr) TBR(TBR(TBR(TBR(tr)))), 7L)
    identical(tr, median(llis, Distance)) &
    identical(tr, median(nnis, Distance)) &
    identical(tr, median(sprs, Distance)) &
    identical(tr, median(tbrs, Distance))
  }, logical(1)))
}

Calc <- function (leaves) {
  vapply(TDFunctions, Successes, integer(1),
         trees = bullseyeTrees[[leaves]][1:5])
}

medianResults <- vapply(paste(c(10, 20), 'leaves'), Calc,
                        integer(length(TDFunctions)))

# Success comes too easily for this to be interesting.
#usethis::use_data(medianResults, overwrite = TRUE)
