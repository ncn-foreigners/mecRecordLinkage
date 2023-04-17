#' @title Supervised Maximum Entropy Classifier for record linkage
#' @author Łukasz Chrostowski
#' @description \code{mecSup} fits MEC record linkage for a supervised learning.
#' @param A The first data set to be linked.
#' @param B The second data set to be linked.
#' @param Omega Training set containing pairs of comparison vectors.
#' @param vars `character` vector of variables that should be compared.
#' @param g The dummy variable for matched and non-matched pairs in the Omega set.
#' @param blockvars The variables defining the blocks or strata for which all pairs of x and y will be generated.
#' @param string_comparator named vector of functions for comparing pairs of records.
#' @param control A list indicating control parameters to use in fitting model.
#' @param prob_ratio Method for the probability ratio estimation.
#'
#'
#' @references D. Lee, L. C. Zhang and J. K. Kim. Maximum entropy classification for record linkage (2022)
#'
#' @importFrom reclin2 compare_pairs
#' @importFrom reclin2 compare_vars
#' @importFrom reclin2 pair_blocking
#' @importFrom reclin2 pair
#' @importFrom FixedPoint FixedPoint
#' @importFrom reclin2 identical
#' @export

mecSup <- function(A,
                   B,
                   Omega = NULL,
                   vars,
                   g,
                   blockvars = NULL,
                   string_comparator = list(reclin2::identical()),
                   control = control_mec(),
                   prob_ratio = "1") {

  nA <- nrow(A)
  nB <- nrow(B)
  fixed_method <- control$fixed_method
  all <- control$all

  if(is.null(blockvars)) {
    pairs <- reclin2::pair(A, B, deduplication = control$deduplication, add_xy = control$add_xy)
    reclin2::compare_pairs(pairs, on = vars, comparators = string_comparator, inplace = TRUE)
  } else {
    pairs <- reclin2::pair_blocking(A, B, on = blockvars, deduplication = control$deduplication, add_xy = control$add_xy)
    reclin2::compare_pairs(pairs, on = vars, comparators = string_comparator, inplace = TRUE)
  }
  pairs[is.na(pairs)] <- FALSE
  gamma <- reclin2::tabulate_patterns(pairs)
  n_pairs <- nrow(pairs)
  K <- length(vars)

  #supervised
  u_loc <- which(g == 0)
  m_loc <- which(g == 1)
  M <- Omega[m_loc, ]
  U <- Omega[u_loc, ]
  nM <- nrow(M)
  n <- nrow(Omega)

  # train a classifier

  pi <- nM/n
  theta <- params_formula(M, nM)
  etaq <- params_formula(Omega, n)
  etau <- params_formula(U, nrow(U))
  m_gamma_omega <- gamma_formula(theta, Omega)
  u_gamma_omega <- gamma_formula(etau, Omega)
  q_gamma_omega <- gamma_formula(etaq, Omega)

  eta <- switch(prob_ratio,
                "1" = etau,
                "2" = etaq)

  r_train <- switch(prob_ratio,
                    "1" = m_gamma_omega/u_gamma_omega,
                    "2" = m_gamma_omega/q_gamma_omega)
  Omega$r <- r_train
  Omega <- Omega[order(-Omega$r), ]

  # classifying
  u_gamma <- gamma_formula(etau, subset(gamma, select = c(-n)))
  m_gamma <- gamma_formula(theta, subset(gamma, select = c(-n)))
  q_gamma <- gamma_formula(etaq, subset(gamma, select = c(-n)))

  r_class <- switch(prob_ratio,
                    "1" = m_gamma/u_gamma,
                    "2" = m_gamma/q_gamma)

  # Application to the target pairs
  #nM_sup <- sum(pi*r_class/(pi * (r_class - 1) + 1) * gamma$n)
  nM_start <- sum(ifelse(apply(pairs, 1, sum) == K, 1, 0))
  fun <- fixed_nM(n_gamma = gamma$n, r_gamma = r_class, n = n_pairs)
  #nM_sup <- spuRs::fixedpoint(ftn = fun, x0 = nM_start)
  nM_sup <- FixedPoint::FixedPoint(Function = fun, Inputs = nM_start, Method = fixed_method)$FixedPoint

  m <- gamma_formula(theta, subset(pairs, select = vars))
  u <- gamma_formula(etau, subset(pairs, select = vars))
  q <- gamma_formula(etaq, subset(pairs, select = vars))
  pairs$r <- switch(prob_ratio,
                    "1" = m/u,
                    "2" = m/q)
  #m/(pi * m + (1 - pi)*u)

  pairs_ord <- pairs[order(-pairs$r), ]
  M <- pairs_ord[1:round(nM_sup), ]
  U <- pairs_ord[-c(1:round(nM_sup)), ]


  max_class_entropy <- class_entropy(round(nM_sup), M$r)
  M$selected <- TRUE
  U$selected <- FALSE
  to_link <- rbind(M, U)
  linked_data <- reclin2::link(subset(to_link, select = c(-r)), selection = "selected", all = all)
  matching_prob <- nM_sup/n_pairs
  params <- t(data.frame(m = theta, u = eta))

  list(params = params,
       matched_ratio = matching_prob,
       max_class_entropy = max_class_entropy,
       M = M,
       U = U,
       linked_data = linked_data)
}
