##' Peak prevalence
##'
##' Approximation of peak prevalence in the SIR model with vital dynamics.
##'
##' The peak prevalence is
##'
##' \deqn{
##'   y_{\rm max} \approx 1-
##'     \dfrac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big)
##'     - \varepsilon 
##'   \Big(1-\dfrac{1}{{\cal R}_0}\Big)
##'   \dfrac{(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) - 1}
##'   {(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) + {\cal R}_0} \,.
##' }
##'
##' The expression is exact for \eqn{\varepsilon = 0}.
##'
##' @inheritParams fizzle_prob
##' @param epsilon mean infectious period as a proportion of mean
##'     lifetime (\eqn{\varepsilon = 0} corresponds to infinite
##'     lifetime, i.e., no mortality)
##'
##' @return real number between 0 and 1
##' @export
##'
##' @examples
##' peak_prev(R0=2)
##' peak_prev(R0=2, epsilon=c(0, 0.001, 0.01, 0.1, 0.9))
##' peak_prev(R0=c(2,4), epsilon=0.1)
##' peak_prev(R0=c(2,4), epsilon=c(0,0.1)) # FIX: fails when both are vectors
##'
peak_prev <- function( R0, epsilon=0 ) {
    pc <- 1 - 1/R0 # p_crit
    lnR <- log(R0)
    ymax <- 1 - (1/R0)*(1 + lnR) - epsilon*pc*(lnR/pc - 1)/(lnR/pc + R0)
    return(ymax)
}
