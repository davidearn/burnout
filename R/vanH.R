##' Persistence probability via van Herwaarden (1997) approximation
##'
##' @seealso \code{\link{burnout_prob_vanH}}, \code{\link{P1_prob_MS}},
##'     \code{\link{P1_prob}}
##' 
##' @inheritParams P1_prob_other
##' @inheritParams stats::integrate
##'
##' @importFrom emdbook lambertW
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Ballard2016}{burnout}
##' 
##' \insertRef{vanH1997}{burnout}
##'
##' @export
##'
P1_prob_vanH <- function(R0, epsilon, k=1, N=10^6, subdivisions=1000L,
                         ... ) {

    P1 <- P1_prob_other(R0 = R0, epsilon = epsilon,
                        burnout_prob_fun = burnout_prob_vanH, k = k, N = N,
                        subdivisions = subdivisions, ... )
    return(P1)
}

##' Burnout probability based on van Herwaarden (1997)
##'
##' @details
##' See equation (8) of Ballard et al (2016)
##'
##' @seealso \code{\link{P1_prob_vanH}}
##'
##' @inheritParams P1_prob_vanH
##' @param tiny amount which to avoid integration limits (where
##'     integrand blow up)
##' 
##' @importFrom emdbook lambertW
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Ballard2016}{burnout}
##' 
##' \insertRef{vanH1997}{burnout}
##' 
##' @export
##'
burnout_prob_vanH <- function( R0, epsilon, N=10^6,
                              subdivisions=1000L,
                              tiny=1e-12,
                              ... ) {
    ## choose units such that gamma+mu = 1, i.e., mean time infected
    ## period is 1:
    beta <- R0
    mu <- epsilon
    gamma <- 1-epsilon

    W0 <- emdbook::lambertW

    bog <- beta/gamma
    x1A <- (-1/bog)*W0(-bog*exp(-bog))

    integrand <- function(s) {
        (x1A/(1-x1A)) * gamma*(s - s*log(s) - 1) /
            ((beta*s^2)*(1-s+(1/bog)*log(s))) +
            1/(s-x1A)
    }
    messy.integral <-
        try(stats::integrate(f=integrand, lower=x1A+tiny, upper=1-tiny,
                      subdivisions=subdivisions, ...)$value)
    ## FIX: this destroys the automatic vectorization:
    if ("try-error" %in% class(messy.integral)) return(NA)

    C3 <- -log(-beta*x1A / (beta*x1A - gamma)) - messy.integral

    K <- (1/mu)*exp(((1/mu)*(beta*x1A + (beta-gamma)*log(1-x1A))) + C3)

    Gamma <- base::gamma

    p0 <- exp(-K*N*mu^2*(beta/mu)^((beta-gamma-mu)/mu) * exp(-beta/mu) /
              ((gamma+mu)*Gamma((beta-gamma-mu)/mu))
              )

    return(p0)
}
