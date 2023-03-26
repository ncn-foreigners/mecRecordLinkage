#' Title
#'
#' @param comp_rate - a
#' @param deduplication - a
#' @param add_xy - a
#' @param maxit - a
#' @param eps - a
#' @param theta_est - a
#' @param target_flr - a
#' @param treshold - a
#' @param theta_start - a
#'
#' @export
control_mec <- function(comp_rate = .9,
                        deduplication = FALSE,
                        add_xy = TRUE,
                        maxit = 100,
                        eps = 1e-3,
                        theta_est = "1",
                        target_flr = .05,
                        treshold = 1, # to consider
                        increase_rate = 4,
                        theta_start = .5){


  list(comp_rate = comp_rate,
       deduplication = deduplication,
       add_xy = add_xy,
       maxit = maxit,
       eps = eps,
       theta_est = theta_est,
       target_flr = target_flr,
       treshold = treshold,
       increase_rate = increase_rate,
       theta_start = theta_start)
}
