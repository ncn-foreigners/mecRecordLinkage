#' Title mec
#'
#' @param A - a
#' @param B - a
#' @param vars - a
#' @param g - a
#' @param blockvars - a
#' @param control - a
#'
#' @importFrom reclin2 compare_pairs
#' @importFrom reclin2 compare_vars
#' @importFrom reclin2 pair_blocking
#' @importFrom reclin2 pair
#' @export
#'
mec <- function(A, B, Omega = NULL, vars, g, blockvars = NULL, control = control_mec()) {

  nA <- nrow(A)
  nB <- nrow(B)

  if(is.null(blockvars)) {
    pairs <- reclin2::pair(A, B, deduplication = control$deduplication, add_xy = control$add_xy)
  } else {
    pairs <- reclin2::pair_blocking(A, B, on = blockvars, deduplication = control$deduplication, add_xy = control$add_xy)
  }

  gamma <- reclin2::tabulate_patterns(pairs)

  #supervised
  u_loc <- which(g == 0)
  m_loc <- which(g == 1)
  M <- Omega[m_loc, ]
  U <- Omega[u_loc, ]
  nM <- nrow(M)
  n <- nrow(Omega)

  ### probability ratio 2

  pi <- nM/n
  theta <- params_formula(M, nM)
  eta <- params_formula(Omega, n)
  m_gamma_omega <- gamma_formula(theta, Omega)
  u_gamma_omega <- gamma_formula(eta, Omega)

  u_gamma <- gamma_formula(eta, subset(gamma, select = c(-n)))
  m_gamma <- gamma_formula(theta_start, subset(gamma, select = c(-n)))
  r_gamma <- m_gamma/u_gamma
  nM <- round(pi*r_gamma/(pi * (r_gamma - 1) + 1) * gamma$n)

  #### probability ratio 1
  m_gamma <- r_gamma/(pi* (r_gamma - 1) + 1)

}
