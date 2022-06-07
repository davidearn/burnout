##' Persistence probability via Meerson and Sasorov (2009) approximation
##'
##' @seealso \code{\link{burnout_prob_MS}}, \code{\link{P1_prob_MS}},
##'     \code{\link{P1_prob}}
##' 
##' @inheritParams P1_prob
##' @inheritParams stats::integrate
##' @param ... additional arguments pass to \code{\link[stats]{integrate}}
##'
##' @importFrom emdbook lambertW
##' @importFrom stats integrate
##'
##' @export
##'
P1_prob_MS <- function( R0, epsilon, k=1, N=10^6, subdivisions=1000L, ... ) {

    ## FIX: this is identical to P1_prob_vanH except that burnout_prob_MS
    ##      replaces burnout_prob_vanH.   This is dumb.

    ## Ballard et al (2016) use the notation p0:
    p0 <- burnout_prob_MS(R0 = R0, epsilon = epsilon, N = N,
                            subdivisions = subdivisions, ...)

    ## FIX: This is identical to the code in P1_prob() except that the
    ##      probability of burning out conditional on not fizzling
    ##      (p0) is calculated above via van H's formulae.
    fizz <- fizzle_prob(R0=R0, k=k)
    ## pk = probability of not fizzling:
    notfizz <- 1 - fizz
    ## probability of not fizzling and then burning out:
    notfizz.and.burn <- notfizz * p0
    ## probability of either fizzling or burning out:
    fizz.or.burn <- fizz + notfizz.and.burn
    ## persist after neither fizzling nor burning out:
    P1 <- 1 - fizz.or.burn
    return(P1)
}

##' Burnout probability based on Meerson and Sasorov (2009)
##'
##' @details
##' See equation (9) of Ballard et al (2016)
##'
##' The derivation assume \eqn{N S_0 \gg 1}, where \eqn{S_0} is
##' defined in equation (9) of Ballard et al (2016).
##'
##' @seealso \code{\link{P1_prob_MS}}
##'
##' @inheritParams P1_prob_MS
##' 
##' @importFrom emdbook lambertW
##' @importFrom stats integrate
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
    W0 <- emdbook::lambertW
    xm <- (-1/R0) * W0((-1/R0)*exp(-1/R0)) - 1

    integrand <- function(s) {
        s*(s+delta) / ((1+s)^2 * (s - (1-delta)*log(1+s))) -
            xm / ((1+xm)*(s-xm))
    }
    Q1 <- try(stats::integrate(f=integrand, lower=0, upper=xm,
                               subdivisions=subdivisions, ...)$value)
    if ("try-error" %in% class(Q1)) return(NA)

    ym <- (delta+xm)*xm/(1+xm)*(-xm/delta)^(K*delta) *
        exp(K*(xm+delta) - (1+1/xm)*Q1)

    C <- ym*delta / (2*pi*(1-delta))

    S0 <- C * sqrt( 2*pi / (K*delta) )

    p0 <- exp(-N * S0)

    return(p0)
}
