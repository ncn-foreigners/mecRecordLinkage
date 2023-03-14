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
mec <- function(A, B, vars, g, blockvars = NULL, string_comparator = FALSE, control = control_mec()) {

  maxit <- control$maxit
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


  subpairs <- subset(pairs, select = vars)
  n <- nrow(subpairs)
  K <- length(vars)
  pi <- nM/nOmega

  gamma <- unique(subpairs)
  gamma <- sa.data.frame(gamma)
  subpairs <- as.data.frame(subpairs)
  theta_start <- rep(1, K)
  M_start <- subpairs[apply(subpairs, 1, all)]
  nM_start <- sum(apply(subpairs, 1, all))
  g <- vector(mode = "numeric", length = length(gamma))

  it <- 1
  while(TRUE) {

    it <- it + 1
    g <- ifelse(apply(M_start, 1, function(x) any(identical(x, subpairs))), 1, 0) # to fix
    eta <- theta_formula(pairs, n)
    u_gamma <- m_gamma(eta)
    theta <- 1/nM * g * apply(pairsM, 2, sum)
    m_gamma <- m_gamma(theta_start)
    r_gamma <- m_gamma/u_gamma
    nM <- sum(g)
    nM <- nMformula(gamma = gamma, r_gamma = r_gamma, n = n, nM = nM)

  }

}


error_rate <- function(nM, gM) {

  flr <- 1/nM * sum(1 - gM)
  mmr <- (1 - 1/nm) * sum(gM)

  list(flr = flr,
       mmr = mmr)
}

nMformula <- function(gamma, pi = NULL, r_gamma, n, nM) { #to consider

  n_gamma <- nrow(gamma)
  g_gamma <- min(nM*r_gamma/(nM * (r_gamma - 1) + n), 1)
  sum(n_gamma*g_gamma)
}

theta_formula <- function(pairsM, nM) {

  theta <- 1/nM * apply(pairsM, 2, sum)
  theta
}

class_entropy <- function(nM, r_gamma) {
  1/nM * sum(log(r_gamma))
}

Obj_function <- function(u_gamma, r_gamms, nMgamma, nM) {
  Q <- sum(u_gamma/nMgamma * r_gamma) - 1/nM * sum(log(r_gamma))
  Q
}

m_gamma <- function(theta){
  theta*(1 - theta)
}


