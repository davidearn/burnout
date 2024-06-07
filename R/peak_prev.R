##' Peak prevalence
##'
##' Approximation of peak prevalence in the SIR model with vital dynamics.
##'
##' The peak prevalence is
##'
##' \deqn{
##'   y_{\rm max} \approx 1-
##'     \frac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big)
##'     - \varepsilon
##'   \Big(1-\frac{1}{{\cal R}_0}\Big)
##'   \frac{(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) - 1}
##'   {(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) + {\cal R}_0} \,.
##' }
##'
##' The expression is exact for \eqn{\varepsilon = 0}.
##'
##' @seealso \code{\link{peak_prev_nvd}}
##'
##' @inheritParams fizzle_prob
##' @param epsilon mean infectious period (\eqn{1/(\gamma+\mu)}) as a proportion of mean
##'     lifetime (\eqn{1/\mu}; \eqn{\varepsilon = 0} corresponds to infinite
##'     lifetime, i.e., no mortality)
##'
##' @return real number between 0 and 1
##' @export
##'
##' @examples
##' peak_prev(R0=2)
##' peak_prev(R0=2, epsilon=c(0, 0.001, 0.01, 0.1, 0.9))
##' peak_prev(R0=c(2,4), epsilon=0.1)
##' peak_prev(R0=2, epsilon=0.001, eta=c(0.01, 0.1))
##' peak_prev(R0=c(2,4), epsilon=c(0,0.1)) # FIX: fails when both are vectors
##'
peak_prev <- function( R0, epsilon=0, eta=0 ) {
    pc <- 1 - 1/R0 # p_crit
    lnR <- log(R0)
    ymax <- 1 - (1/R0)*(1 + lnR) -
       (epsilon+eta)*pc*(lnR/pc - 1)/(lnR/pc + R0) +
       (eta/R0)*(lnR - pc)
    return(ymax)
}

##' Peak prevalence
##'
##' Exact peak prevalence in the SIR model with \emph{no vital dynamics} (nvd).
##'
##' \deqn{
##'   y_{\rm max} = 1 -
##'     \frac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big) \,.
##' }
##'
##' The optional argument (\code{epsilon}) is ignored if given.  The
##' reason for including it is so this function can be used as an
##' alternative to the more accurate \code{\link{peak_prev}} function
##' that does depend on \eqn{\varepsilon}.
##'
##' @seealso \code{\link{peak_prev}}, \code{\link{x_in}}
##'
##' @inheritParams peak_prev
##'
##' @return real number between 0 and 1
##' @export
##'
peak_prev_nvd <- function( R0, epsilon=0 ) {
    ymax <- 1 - (1/R0)*(1 + log(R0))
    return(ymax)
}
