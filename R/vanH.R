##' Persistence probability via van Herwaarden (1997) approximation
##'
##' @seealso \code{\link{burnout_prob_vanH}}, \code{\link{P1_prob}}
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
P1_prob_vanH <- function( R0, epsilon, k=1, N=10^6, subdivisions=1000L, ... ) {

    ## Ballard et al (2016) use the notation p0:
    p0 <- burnout_prob_vanH(R0 = R0, epsilon = epsilon, N = N,
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

##' Burnout probability based on van Herwaarden (1997)
##'
##' @details
##' See equation (8) of Ballard et al (2016)
##'
##' @seealso \code{\link{P1_prob_vanH}}
##'
##' @inheritParams P1_prob_vanH
##' 
##' @importFrom emdbook lambertW
##' @importFrom stats integrate
##'
##' @export
##'
burnout_prob_vanH <- function( R0, epsilon, N=10^6, subdivisions=1000L, ... ) {
    ## choose units such that gamma+mu = 1, i.e., mean time infected is 1:
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
        try(stats::integrate(f=integrand, lower=x1A, upper=1,
                      subdivisions=subdivisions, ...)$value)
    if ("try-error" %in% class(messy.integral)) return(NA)

    C3 <- -log(-beta*x1A / (beta*x1A - gamma)) - messy.integral

    K <- (1/mu)*exp(((1/mu)*(beta*x1A + (beta-gamma)*log(1-x1A))) + C3)

    Gamma <- base::gamma

    p0 <- exp(-K*N*mu^2*(beta/mu)^((beta-gamma-mu)/mu) * exp(-beta/mu) /
              ((gamma+mu)*Gamma((beta-gamma-mu)/mu))
              )

    return(p0)
}
