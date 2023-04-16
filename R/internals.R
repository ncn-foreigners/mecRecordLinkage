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

fixed_nM <- function(n_gamma, n, r_gamma){
  function(nM){
    sum(n_gamma * nM*r_gamma/(nM*(r_gamma - 1) + n))
  }
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

mlee <- function(nM, n, A, B, M, u, Omega){

  if (nA > nB) {
    C <- B
    B <- A
    A <- B
    nA <- nrow(A)
    nB <- nrow(B)
  }
  #zA <- apply(A, 1, unique)
  #zB <- apply(B, 1, unique)

  pi <- nM/n
  p <- nA * pi
  mA <- apply(A, 2, table)
  nD <- length(unlist(mA))
  mle1 <- function(A){
    mA <- apply(A, 2, table)
    function(x){
      idx <- vector(mode = "numeric", length = length(mA))
      logs <- vector(mode = "numeric", length = length(mA))
      start <- 1
      a <- 0
      for (i in 1:length(mA)) {
        idx[i] <- length(mA[[i]]) + a
        logs[i] <- log(sum(x[start:idx[i]] * mA[[i]]))
        start <- idx[i] + 1
        a <- idx[i]
      }
      sum(logs)
    }
  }


  mle2 <- function(B, p, deltaB){ # EM algorithm here
    mB <- apply(B, 2, table)
    uB <- apply(B, 2, table)
    function(x){
      #log_likeB <- sum(delta * log(p) * prod(x * mB)) + sum((1 - delta) * log(1 - p) * prod(x * uB)) # x is mkd and ukd
      sum(delta * unlist(lapply(mA, function(z) log(p*prod(mapply("*", x, z)))))) + sum((1 - delta) * unlist(lapply(mA, function(z) log((1 - p)*prod(mapply("*", x, z))))))
      for (i in 1:length(mB)) {
        idx[i] <- length(mA[[i]]) + a
        mbk <- sum(x[start:idx[i]] * mB[[i]])
        ubk <- sum(x[start:idx[i]] * uB[[i]])
        start <- idx[i] + 1
        a <- idx[i]
      }
      deltaB * log(p * prod(mbk)) + (1 - deltaB) * log((1 - p) * prod(ubk))
    }
  }
  #eta <- ((1 - p) * sum(u*m) + p * (1 - 1/nA) * sum(m^2))/(1 - p/nA)
  #eta <- apply(D, 1, function(x) (1 - p) * sum())
  #eta
  m
}
