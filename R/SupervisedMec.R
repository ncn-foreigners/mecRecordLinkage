#' Title mec
#'
#' @param A - a
#' @param B - a
#' @param vars - a
#' @param g - a
#' @param blockvars - a
#' @param string_comparator - a
#' @param control - a
#'
#' @importFrom reclin2 compare_pairs
#' @importFrom reclin2 compare_vars
#' @importFrom reclin2 pair_blocking
#' @importFrom reclin2 pair
#' @export
#'
mec <- function(A, B, Omega = NULL, vars, g, blockvars = NULL, string_comparator = FALSE, control = control_mec()) {

  nA <- nrow(A)
  nB <- nrow(B)

  if(is.null(blockvars)) {
    pairs <- reclin2::pair(A, B, deduplication = control$deduplication, add_xy = control$add_xy)
  } else {
    pairs <- reclin2::pair_blocking(A, B, on = blockvars, deduplication = control$deduplication, add_xy = control$add_xy)
  }
  if (string_comparator) {
    reclin2::compare_pairs(pairs, on = blockvars, inplace = TRUE, default_comparator = jaro_winkler(treshold = control$comp_rate))
  } else {
    reclin2::compare_pairs(pairs, on = blockvars, inplace = TRUE)
  }


  #supervised
  u_loc <- which(g == 0)
  m_loc <- which(g == 1)
  M <- Omega[m_loc, ]
  U <- Omega[u_loc, ]

}
