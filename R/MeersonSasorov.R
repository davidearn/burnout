##' Persistence probability via Meerson and Sasorov (2009) approximation
##'
##' @seealso \code{\link{burnout_prob_MS}}, \code{\link{P1_prob_MS}},
##'     \code{\link{P1_prob}}
##' 
##' @inheritParams P1_prob_other
##'
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Meerson2009}{burnout}
##' 
##' @export
##'
P1_prob_MS <- function( R0, epsilon, k=1, N=10^6, subdivisions=1000L, ... ) {

    P1 <- P1_prob_other(R0 = R0, epsilon = epsilon,
                        burnout_prob_fun = burnout_prob_MS, k = k, N = N,
                        subdivisions = subdivisions, ... )
    return(P1)
}

##' Burnout probability based on Meerson and Sasorov (2009)
##'
##' @details
##' See equation (9) of Ballard et al (2016)
##'
##' The derivation assumes \eqn{N S_0 \gg 1}, where \eqn{S_0} is
##' defined in equation (9) of Ballard et al (2016).
##'
##' @seealso \code{\link{P1_prob_MS}}
##'
##' @inheritParams P1_prob_MS
##' 
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Meerson2009}{burnout}
##' 
##' @export
##'
burnout_prob_MS <- function( R0, epsilon, N=10^6, subdivisions=1000L, ... ) {
    ## choose units such that gamma+mu = 1, i.e., mean time infected is 1:
    beta <- R0
    mu <- epsilon
    gamma <- 1-epsilon

    K <- beta/mu
    delta <- 1 - 1/R0
    xm <- (-1/R0) * W0(-R0*exp(-R0)) - 1

    integrand <- function(s) {
        s*(s+delta) / ((1+s)^2 * (s - (1-delta)*log(1+s))) -
            xm / ((1+xm)*(s-xm))
    }
    Q1 <- try(stats::integrate(f=integrand, lower=0, upper=xm,
                               subdivisions=subdivisions, ...)$value)
    ## FIX: this destroys the automatic vectorization:
    if ("try-error" %in% class(Q1)) return(NA)

    ym <- (delta+xm)*xm/(1+xm)*(-xm/delta)^(K*delta) *
        exp(K*(xm+delta) - (1+1/xm)*Q1)

    C <- ym*delta / (2*pi*(1-delta))

    S0 <- C * sqrt( 2*pi / (K*delta) )

    p0 <- exp(-N * S0)

    return(p0)
}
