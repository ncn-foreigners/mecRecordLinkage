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

params_formula <- function(set, nM) {

  par <- 1/nM * apply(set, 2, sum)
  par
}

params_formula_theta <- function(set, nM, g_gamma) {

  par <- 1/nM * t(as.matrix(g_gamma)) %*% as.matrix(set)
  as.vector(par)
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

mle <- function(){




}
