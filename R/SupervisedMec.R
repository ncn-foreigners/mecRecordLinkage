#' Title mec
#'
#' @param A
#' @param B
#' @param vars
#' @param blockvars
#'
#' @importFrom reclin2 compare_pairs
#' @importFrom reclin2 compare_vars
#' @importFrom reclin2 pair_blocking
#' @importFron reclin2 pair
#' @export
#'
#' @examples
mec <- function(A, B, Omega = NULL, vars, g, blockvars = NULL, string_comparator = FALSE, control = control_mec()) {

  nA <- nrow(A)
  nB <- nrow(B)

  if(is.null(blockvars)) {
    pairs <- pair(A, B, deduplication = control$deduplication, add_xy = control$add_xy)
  } else {
    pairs <- pair_blocking(A, B, on = blockvars, deduplication = control$deduplication, add_xy = control$add_xy)
  }
  if (string_comparator) {
    compare_pairs(pairs, on = blockvars, inplace = TRUE, default_comparator = jaro_winkler(treshold = control$comp_rate))
  } else {
    compare_pairs(pairs, on = blockvars, inplace = TRUE)
  }

  #supervised
  gU <- which(g == 0)
  gM <- which(g == 1)
  gammaA <- pairs[gM,]
  gammaU <- pairs[gU,]

}
