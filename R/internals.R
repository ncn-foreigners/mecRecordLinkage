# Internals function

error_rate <- function(nM, g) {

  flr <- 1/nM * sum(1 - g)
  mmr <- 1 - sum(g/nM)

  list(flr = flr,
       mmr = mmr)
}

nMformula <- function(gamma, pi = NULL, r_gamma, n, nM) {

  n_gamma <- gamma$n
  g_gamma <- sapply(nM*r_gamma/(nM * (r_gamma - 1) + n), function(x) min(c(x, 1)))
  sum(n_gamma*g_gamma)
}

params_formula <- function(set, n) {

  par <- 1/n * apply(set, 2, sum)
  par
}

params_formula_theta <- function(set, nM, g_gamma) {

  #par <- 1/nM * t(as.matrix(g_gamma)) %*% as.matrix(set)
  par <- 1/nM * apply(set, 2, function(x) sum(g_gamma * x))
  names(par) <- names(set)
  par
}

#params_formula_eta <- function(pi, nA, nB, )

class_entropy <- function(nM, r_gamma) {
  1/nM * sum(log(r_gamma))
}

Obj_function <- function(u_gamma, r_gamms, nMgamma, nM) {
  Q <- sum(u_gamma/nMgamma * r_gamma) - 1/nM * sum(log(r_gamma))
  Q
}

gamma_formula <- function(par, subpairs){
  apply(subpairs, 1, function(x) prod(par^x * (1 - par)^(1-x)))
}

mle <- function(nM, n, nA, nB, A, B, M, Omega){

  if (nA > nB) {
    C <- B
    B <- A
    A <- B
    nA <- nrow(A)
    nB <- nrow(B)
  }

  pi <- nM/n
  p <- nA * pi
  e <- rbinom(n = nM, size = 1, prob = 1/2)

  mA <- unlist(apply(A, 2, table))
  log_likeA <- function(x) sum(log(x * mA)) # x is mkd

  mB <- unlist(apply(B, 2, table))
  uB <- unlist(apply(B, 2, table))
  log_likeB <- function(x) sum(delta * log(p) * prod(x * mB)) + sum((1 - delta) * log(1 - p) * prod(x * uB)) # x is mkd and ukd

  um <- u*m

  eta <- ((1 - p) * sum(um) + p * (1 - 1/nA) * sum(m^2))/(1 - p/nA)

  eta <- aplly(D, 1, function(x) (1 - p) * sum())



}
