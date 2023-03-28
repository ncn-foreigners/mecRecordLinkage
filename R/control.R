#' @title Control parameters Maximum entropy classification
#' @author Łukasz Chrostowski
#'
#' #' @description \code{control_mec} constructs a list with all necessary control parameters
#' for MEC record linkage.
#'
#' @param deduplication generate pairs from only x. Ignore y. This is useful for deduplication of x. Default is `FALSE`.
#' @param add_xy add x and y as attributes to the returned pairs. This makes calling some subsequent operations that need x and y. Default is `FALSE`.
#' @param maxit The number of maximum iteration for fitting model in the unsupervised learning. Default is `100`.
#' @param eps Tolerance for fitting algorithms by default \code{1e-3}.
#' @param theta_est formula for theta estimation. Default is `1`.
#' @param target_flr The target error rate (false link rate) for the fitting. Default is `.05`.
#' @param treshold The threshold value for fitting the algorithm by the error rates. Default is `1`.
#' @param increase_rate The increasing value for fitting the algorithm by the error rates. Default is `4`.
#' @param theta_start starting value for the theta parametet. Default is `.5`.
#'
#' @export
control_mec <- function(deduplication = FALSE,
                        add_xy = TRUE,
                        maxit = 100,
                        eps = 1e-3,
                        theta_est = "1",
                        target_flr = .05,
                        treshold = 1, # to consider
                        increase_rate = 4,
                        theta_start = .5){


  list(deduplication = deduplication,
       add_xy = add_xy,
       maxit = maxit,
       eps = eps,
       theta_est = theta_est,
       target_flr = target_flr,
       treshold = treshold,
       increase_rate = increase_rate,
       theta_start = theta_start)
}
