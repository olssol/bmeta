#' Call STAN models
#'
#' Call STAN models. Called by \code{psrwe_powerp}.
#'
#' @param lst_data List of study data to be passed to STAN
#' @param stan_mdl STAN model name
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#' @return Result from STAN sampling
#'
#' @export
#'
bmeta_stan <- function(lst_data,
                       stan_mdl = c("riskdiff", "riskratio"),
                       chains = 4, iter = 2000, warmup = 1000,
                       control = list(adapt_delta = 0.95), ...) {

    stan_mdl <- match.arg(stan_mdl)
    stan_rst <- rstan::sampling(stanmodels[[stan_mdl]],
                                data    = lst_data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...)

    stan_rst
}

#' Risk difference or Risk Ratio
#'
#'
#' @export
#'
bmeta_ana <- function(ctl_event, ctl_total, trt_event, trt_total,
                      type = c("riskdiff", "riskratio"),
                      tau_max = 2, tau2_eta = 100,
                      ...) {

    type     <- match.arg(type)
    ns       <- length(ctl_event)
    lst_data <- list(ns       = ns,
                     n_pla    = cbind(ctl_event, ctl_total),
                     n_trt    = cbind(trt_event, trt_total),
                     tau_max  = tau_max,
                     tau2_eta = tau2_eta)

    rst_stan <- bmeta_stan(lst_data, stan_mdl = type, ...)

    rst <- switch(type,
                  riskdiff  = bmeta_diff(rst_stan),
                  riskratio = bmeta_ratio(rst_stan))

    rst        <- c(rst, list(rst_stan = rst_stan))
    class(rst) <- "bmeta"

    invisible(rst)
}

#' Print
#'
#' @method print bmeta
#'
#' @export
#'
print.bmeta <- function(x, ...) {
    cat("Type:     ", x$type, "\n")
    cat("Estimate: ",
        x$Estimate[1],
        "(", x$Estimate[2], ",",
        x$Estimate[3], ")",
        sep = "")

    cat("\n")
    cat("M1      : ", x$M1, "\n")
    cat("M2      : ", x$M2, "\n")
}

#'
#'
#'
#' @export
#'
bmeta_ratio <- function(rst_stan) {

    eta <- rstan::extract(rst_stan, pars = "eta")$eta
    eta <- exp(eta)
    m   <- mean(eta)
    lb  <- quantile(eta, 0.025)
    ub  <- quantile(eta, 0.975)

    m1  <- 1 / ub
    m2  <- exp(0.5 * log(m1))

    list(eta        = eta,
         type       = "Risk Ratio",
         Estimate   = c(estimate = m, lb = lb, ub = ub),
         M1         = m1,
         M2         = m2)
}

#'
#'
#'
#' @export
#'
bmeta_diff <- function(rst_stan) {

    eta <- rstan::extract(rst_stan, pars = "eta")$eta
    m   <- mean(eta)
    lb  <- quantile(eta, 0.025)
    ub  <- quantile(eta, 0.975)

    m1  <- - ub
    m2  <- 0.5 * m1

    list(eta        = eta,
         type       = "Risk Difference",
         Estimate   = c(estimate = m, lb = lb, ub = ub),
         M1         = m1,
         M2         = m2)
}
