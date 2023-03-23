control_mec <- function(comp_rate = .9,
                        deduplication = FALSE,
                        add_xy = TRUE,
                        maxit = 100,
                        eps = 1e-6,
                        theta_est = "1",
                        target_flr = .01,
                        treshold = 1,
                        theta_start = .5){


  list(comp_rate = comp_rate,
       deduplication = deduplication,
       add_xy = add_xy,
       maxit = maxit,
       eps = eps,
       theta_est = theta_est,
       target_flr = target_flr,
       treshold = treshold,
       theta_start = theta_start)
}
