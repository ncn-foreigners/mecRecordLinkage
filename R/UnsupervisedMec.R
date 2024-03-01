#' @title Unsupervised Maximum Entropy Classifier for a record linkage
#' @author Łukasz Chrostowski
#' @description \code{mec} fits MEC record linkage for an unsupervised learning.
#' @param A first data set to be linked.
#' @param B second data set to be linked.
#' @param vars character vector of variables that should be compared.
#' @param blockvars the variables defining the blocks or strata for which all pairs of x and y will be generated.
#' @param error_rate if `TRUE`, estimation is guided by errors rate is generated.
#' @param string_comparator named vector of functions for comparing pairs of records.
#' @param control A list indicating control parameters to use in fitting model.
#'
#' @references D. Lee, L. C. Zhang and J. K. Kim. Maximum entropy classification for record linkage (2022)
#'
#' @importFrom reclin2 compare_pairs
#' @importFrom reclin2 compare_vars
#' @importFrom reclin2 pair_blocking
#' @importFrom reclin2 pair
#' @importFrom reclin2 tabulate_patterns
#' @importFrom reclin2 link
#' @importFrom reclin2 identical
#' @export

mec <- function(A,
                B,
                vars,
                blockvars = NULL,
                error_rate = FALSE,
                string_comparator = list(reclin2::identical()),
                control = control_mec()) {

  maxit <- control$maxit
  eps <- control$eps
  theta_est <- control$theta_est
  eta_est <- control$eta_est
  treshold <- control$treshold
  target_flr <- control$target_flr
  increase_rate <- control$increase_rate
  all <- control$all
  nA <- nrow(A)
  nB <- nrow(B)

  if(is.null(blockvars)) {
    pairs <- reclin2::pair(A, B, deduplication = control$deduplication, add_xy = control$add_xy)
    reclin2::compare_pairs(pairs, on = vars, comparators = string_comparator, inplace = TRUE)
    } else {
    pairs <- reclin2::pair_blocking(A, B, on = blockvars, deduplication = control$deduplication, add_xy = control$add_xy)
    reclin2::compare_pairs(pairs, on = vars, comparators = string_comparator, inplace = TRUE)
    }

  pairs[is.na(pairs)] <- FALSE
  Omega <- subset(pairs, select = vars)
  gamma <- reclin2::tabulate_patterns(pairs)
  n <- nrow(pairs)
  K <- length(vars)
  theta_start <- rep(control$theta_start, K)


  g_start <- ifelse(apply(Omega, 1, sum) == K, 1, 0)

  nM_start <- sum(g_start)
  u_loc = which(g_start==0)
  m_loc = which(g_start==1)

  eta <- params_formula(Omega, n)
  u_gamma <- gamma_formula(eta, subset(gamma, select = c(-n)))
  u_gamma_omega <- gamma_formula(eta, Omega)
  M_start <- Omega[m_loc, ]

  if (!error_rate) {
    if (nM_start == 0) {
      stop("There is no records with all agreements on the key variables. Please provide relevant datasets.")
    }
    it <- 0
    while(TRUE) {
      it <- it + 1
      m_gamma <- gamma_formula(theta_start, subset(gamma, select = c(-n)))
      if (eta_est == "2") {
        u_gamma <- gamma_formula(eta, subset(gamma, select = c(-n)))
      }
      r_gamma <- m_gamma/u_gamma

      nM_prev <- nM_start
      nM <- nMformula(gamma = gamma, r_gamma = r_gamma, n = n, nM = nM_start)

      if (eta_est == "2") {
        u_gamma_omega <- gamma_formula(eta, Omega)
      }
      pairs$r <- gamma_formula(theta_start, Omega)/u_gamma_omega
      pairs_ord <- pairs[order(-pairs$r), ]
      M <- pairs_ord[1:nM, ]
      U <- pairs_ord[-c(1:nM), ]

      M_start <- subset(M, select = vars)
      nM_start <- nM

      if (theta_est == "1") {
        theta <- params_formula(M_start, nM_start)
      } else if (theta_est == "2") {
        g_gamma <- sapply(nM_prev*pairs$r/(nM_prev * (pairs$r - 1) + n), function(x) min(x, 1)) # problem here
        theta <- params_formula_theta(Omega, nM_start, g_gamma)
      }

     if (eta_est == "2") {
        eta <- mle(nM = nM_start, n = n, A = A, B = B, M = M, u, Omega = Omega)
      }

      ##############
      if (sqrt(sum((theta - theta_start)^2)) < eps) break;
      if (round(nM_start) == round(nM_prev)) break;
      theta_start <- theta
    }
  } else {
    nM_prev <- 0
    it <- 0
    while(TRUE) {
      it <- it + 1

      r_gamma <- gamma_formula(theta_start, Omega)/u_gamma_omega
      pairs$r <- r_gamma
      pairs_ord <- pairs[order(-pairs$r),]
      M <- subset(pairs_ord, r >= treshold)
      nM <- nrow(M)
      U <- pairs_ord[-c(1:nM), ]

      r_gamma2 <- gamma_formula(theta_start, subset(gamma, select = c(-n)))/u_gamma
      nM2 <- nMformula(gamma = gamma, r_gamma = r_gamma2, n = n, nM = nM)

      g_gamma <- sapply(nM*M$r/(nM * (M$r - 1) + n), function(x) min(c(x, 1)))
      flr <- error_rate(nM, g_gamma)$flr
      mmr <- error_rate(nM2, g_gamma)$mmr

      M_start <- subset(M, select = vars)

      if (theta_est == "1") {
        theta <- params_formula(M_start, nM)
      } else if (theta_est == "2") {
        g_gamma <- sapply(nM*pairs$r/(nM * (pairs$r - 1) + n), function(x) min(x, 1))
        theta <- params_formula_theta(Omega, nM2, g_gamma)
      }

      if (flr > target_flr) {
        treshold <- treshold + increase_rate # to consider
      } else {
        treshold <- treshold - increase_rate
      }

      ##############
      if (sqrt(sum((theta - theta_start)^2)) < eps) break;
      if (round(nM2) == round(nM_prev)) break;
      #if (abs(flr - target_flr) < eps) break;

      theta_start <- theta
      nM_prev <- nM2
    }
  }

  matching_prob <- nM/n
  max_class_entropy <- class_entropy(nM, M$r)
  M$selected <- TRUE
  U$selected <- FALSE
  to_link <- rbind(M, U)
  linked_data <- reclin2::link(subset(to_link, select = c(-r)), selection = "selected", all = all)
  params <- t(data.frame(m = theta, u = eta))

  structure(
    list(params = params,
       matched_ratio = matching_prob,
       max_class_entropy = max_class_entropy,
       M = M,
       U = U,
       linked_data = linked_data),
    class = "MecRecordLinkage")
}
