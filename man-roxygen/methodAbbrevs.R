#' @details
#' # Methods tested
#'
#' - `dpi`: Different Phylogenetic Information (Smith, forthcoming)
#' - `msid`: Matching Split Information Distance (Smith, forthcoming)
#' - `cid`: Clustering Information Distance (Smith, forthcoming)
#' - `qd`: Quartet divergence (Smith, 2019)
#' - `nts`: Nye _et al._ tree similarity (Nye _et al._ 2006)
#' - `ja2`, `ja4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'            conflicting pairings prohobited
#' - `jna2`, `jna4`: Jaccard-Robinson-Foulds distances with _k_ = 2, 4,
#'            conflicting pairings permitted
#' - `msd`: Matching Split Distance (Bogdanowicz & Giaro, 2012)
#' - `mast`: Size of Maximum Agreement Subtree (Valiente 2009)
#' - `masti`: Information content of Maximum Agreement Subtree
#' - `nni_l`, `nni_t`, `nni_u`: Lower bound, tight upper bound, and upper bound
#'         for  nearest-neighbour interchange distance (Li _et al._ 1996)
#' - `spr`: Approximate SPR distance
#' - `tbr_l`, `tbr_u`: Lower and upper bound for tree bisection and reconnection
#'          (TBR) distance
#' - `rf`: Robinson-Foulds distance (Robinson & Foulds 1985)
#' - `rfi`: Information-corrected Robinson-Foulds distance (Smith, forthcoming)
#' - `path`: Path distance (Steel & Penny 1993), unnormalized
#'
#' @references
#' \insertRef{Bogdanowicz2012}{TreeDist}
#'
#' \insertRef{Li1996}{TreeDist}
#'
#' \insertRef{Kendall2016}{TreeDistData}
#'
#' \insertRef{Nye2006}{TreeDist}
#'
#' \insertRef{Robinson1981}{TreeDist}
#'
#' \insertRef{SmithTern}{TreeSearch}
#'
#' \insertRef{SmithDist}{TreeDist}
#'
#' \insertRef{Steel1993}{TreeDist}
#'
#' \insertRef{Valiente2009}{TreeDist}
