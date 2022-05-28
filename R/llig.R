##' log lower incomplete gamma function
##'
##' Compute the logarithm of the lower incomplete gamma function via
##' \code{\link[stats]{pgamma}}.
##'
##' The \emph{lower incomplete gamma function} ([NIST
##' 8.2.1](https://dlmf.nist.gov/8.2#E1)) is
##'
##' \deqn{\texttt{lig}(a,z) = \int_0^z t^{a-1}e^{-t}\,{\rm d}t \,,
##' \qquad \Re{a}>0.}
##'
##' We are interested only in real arguments, and we write
##' \eqn{\texttt{lig}(a,x)} rather than the usual \eqn{\gamma(a,x)}
##' for this function to avoid confusion with the recovery rate in the
##' SIR model.
##'
##' Since \eqn{\Gamma(a) = \texttt{lig}(a,\infty)}, the
##' \emph{normalized} area under the curve \eqn{t^{a-1}e^{-t}} is
##' 
##' \deqn{P(a,x) = \frac{\texttt{lig}(a,x)}{\Gamma(a)}}
##' 
##' This is the cumulative distribution of the gamma distribution
##' (\code{\link[stats]{pgamma}}), which we can calculate and express
##' on the log scale via
##' 
##' \deqn{\code{pgamma(shape=a, q=x, lower.tail = FALSE, log.p = TRUE)}.}
##' 
##' Since \eqn{\log(\Gamma(a))} can be calculated directly on the
##' log scale via \code{\link[base]{lgamma}}, we can do the entire
##' calculation on the log scale via
##'
##' \deqn{\texttt{llig}(a,z) = P(a,z) + \log(\Gamma(a))}
##' 
##' This approach allows us to compute \code{llig} for very large
##' values of \code{a} and \code{x} without overflow.
##' 
##' Any time you need to convert to values of \eqn{\texttt{lig}(a,x) >
##' 10^{308}} (i.e., in the non-log scale) you're going to be in
##' trouble, although there are packages like
##' \code{\link[Brobdingnag]{Brobdingnag}} that handle these kinds of
##' large numbers.  For our purposes with this package, we can use
##' logspace addition and subtraction.
##'
##' There are other implementations of the lower incomplete gamma
##' function, which fail for large arguments (e.g.,
##' \code{\link[expint]{gammainc}(a,x)},
##' \code{\link[gsl]{gamma_inc}(a,x)}).
##'
##' @param a shape parameter (non-negative real number)
##' @param x quantile parameter (non-negative real number)
##'
##' @return real number
##' @export
##' 
##' @examples
##' llig( a = 1, x = 1 )
##' llig( a = 10^200, x = 10^200 )
##' 
llig <- function(a, x) {
   pgamma(shape = a, q = x, lower.tail = FALSE, log.p = TRUE) +
     lgamma(a)
}
