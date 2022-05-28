##' log lower incomplete gamma function
##'
##' Compute the logarithm of the lower incomplete gamma function via
##' \code{\link[stats]{pgamma}}.
##'
##' The \emph{lower incomplete gamma function} is ([NIST
##' 8.2.1](https://dlmf.nist.gov/8.2#E1)) \deqn{\gamma(a,z) = \int_0^z
##' t^{a-1}e^{-t}\,{\rm d}t \,, \qquad \Re{a}>0.}  We use
##' \eqn{\mathcal{g}} rather than \eqn{\gamma} since \eqn{\gamma} has
##' a standard meaning in the SIR model.
##' 
##' The normalized lower incomplete gamma function is ([NIST
##' 8.2.4](https://dlmf.nist.gov/8.2#E4)) \deqn{P(a,z) =
##' \frac{\gamma(a,z)}{\Gamma(a)}}
##' 
##' \code{\link[stats]{pgamma}} is the \emph{normalized} lower
##' incomplete gamma function.  The normalization constant is given by
##' \code{\link[base]{lgamma}}, which is \eqn{\log(\Gamma(x))}.
##' 
##' \code{pgamma(shape, q, lower.tail = FALSE, log.p = TRUE)} is the
##' \emph{normalized} area under the curve, calculated and expressed
##' on the log scale.  \code{lgamma(a)} is the normalization constant.
##' 
##' This approach allows us to compute \code{lig} for very large
##' values of its arguments (\code{a} and \code{x}) without overflow.
##' 
##' Any time you need to convert to values of \eqn{{\mathcal g} >
##' 10^{308}} (i.e., in the non-log scale) you're going to be in
##' trouble, although there are packages like \code{"Brobdingnag"}
##' that handle these kinds of large numbers.
##'
##' For our present purposes we can use logspace addition and
##' subtraction.
##'
##' @param a shape parameter (non-negative real number)
##' @param x q
##'
##' @return real number
##' @export
##' 
##' @examples
##' llig( a = 1, x = 1 )
##' llig( a = 10^8, x = 10^16 )
##' 
llig <- function(a, x) {
   pgamma(shape = a, q = x, lower.tail = FALSE, log.p = TRUE) +
     lgamma(a)
}
