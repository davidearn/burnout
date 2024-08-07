##' Kendall's \eqn{q} (exact expression)
##'
##' Probability of eventual extinction in a birth-death process
##' starting from one individual.
##'
##' The theorem of \insertCite{Kendall1948b;textual}{burnout} applied
##' to the stochastic SIR model with one infected individual yields
##' 
##' \deqn{
##' 	q({\cal R}_0,\varepsilon,x_{\rm in}) =
##' \left(1+\frac{\varepsilon}{e^{\frac{{\cal R}_{0}}{\varepsilon} (1-x_{\rm in})}
##' \left(\frac{{\cal R}_{0}}{\varepsilon}
##' (1-x_{\rm in})\right)^{-\frac{{\cal R}_{0}}{\varepsilon}\left(1-\frac{1}{{\cal R}_{0}}\right)}
##' {\cal g}\left(\frac{{\cal R}_{0}}{\varepsilon}\left(1-\frac{1}{{\cal R}_{0}}\right),\frac{{\cal R}_{0}}{\varepsilon}(1-x_{\rm in})\right)}\right)^{-1}.
##' }
##'
##'
##' @seealso \code{\link{x_in}}, \code{\link{llig}}, \code{\link{q_approx}}
##'
##' @inheritParams peak_prev
##' @param xin by default set via the approximation coded in \code{\link{x_in}}
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##' 
##' @export
##'
##' @examples
##' q_exact(R0=2, epsilon=0.01)
##' q_exact(R0=2, epsilon=0)
##' q_exact(R0=c(2,4), epsilon=0.01)
##' q_exact(R0=2, epsilon=c(0.01, 0.1))
##' curve(q_exact(x,epsilon=0.001), from=1.01, to=2, las=1, n=1001, ylim=c(0,1))
##' curve(q_exact(x,epsilon=0.01), from=1.01, to=2, las=1, add=TRUE, col="magenta", n=1001)
##' curve(q_exact(x,epsilon=0.1), from=1.01, to=2, las=1, add=TRUE, col="cyan", n=1001)
##'
q_exact <- function(R0, epsilon, eta=0, xin = x_in(R0,epsilon,eta)) {
    eps.plus.eta <- epsilon + eta
    a <- (R0/eps.plus.eta)*(1-1/R0)
    x <- (R0/eps.plus.eta)*(1-xin)
    ##denom <- exp(x) * x^(-a) * lig(a,x)
    log.denom <- x - a*log(x) + llig(a,x)
    denom <- exp(log.denom)
    q <- 1 / (1 + eps.plus.eta/denom)
    return(q)
}

##' Kendall's \eqn{q} (Laplace approximation)
##'
##' Approximate probability of eventual extinction in a birth-death process
##' starting from one individual.
##'
##' Exploiting [Laplace's
##' method](https://en.wikipedia.org/wiki/Laplace%27s_method) we
##' obtain
##' 
##' \deqn{
##' q({\cal R}_0,\varepsilon,x_{\rm in}) \approx
##' \left(1 +
##' \frac{1}{1+\sqrt{\frac{2\pi}{\varepsilon ({\cal R}_{0}-1)}} \;e^{\frac{{\cal R}_{0}}{\varepsilon}\big(\frac{1}{{\cal R}_{0}}-x_{\rm in}\big)}
##' \Big(\frac{1-\frac{1}{{\cal R}_{0}}}{1-x_{\rm in}}\Big)^{\frac{{\cal R}_{0}}{\varepsilon}\big(1-\frac{1}{{\cal R}_{0}}\big)}}
##' \right)^{-1}\,.
##' }
##'
##'
##' @seealso \code{\link{x_in}}, \code{\link{llig}}, \code{\link{q_exact}}
##'
##' @inheritParams q_exact
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##' 
##' @export
##'
##' @examples
##' q_approx(R0=2, epsilon=0.01)
##' q_approx(R0=2, epsilon=0)
##' q_approx(R0=c(2,4), epsilon=0.01)
##' q_approx(R0=2, epsilon=c(0.01, 0.1))
##' curve(q_approx(x,epsilon=0.001), from=1.01, to=2, las=1, n=1001, ylim=c(0,1))
##' curve(q_approx(x,epsilon=0.01), from=1.01, to=2, las=1, add=TRUE, col="magenta", n=1001)
##' curve(q_approx(x,epsilon=0.1), from=1.01, to=2, las=1, add=TRUE, col="cyan", n=1001)
##'
q_approx <- function(R0, epsilon, eta=0, xin = x_in(R0,epsilon,eta)) { 
    eps.plus.eta <- epsilon + eta
    a <- (R0/eps.plus.eta)*(1 - 1/R0)
    b <- (R0/eps.plus.eta)*(1/R0 - xin)
    log.fac1 <- (1/2)*(log(2*pi) - log(eps.plus.eta*(R0-1)))
    log.fac2 <- a*(log(1-1/R0) - log(1-xin))
    log.messy <- log.fac1 + log.fac2 + b
    denom <- 1 + exp(log.messy)
    q <- (1 + 1/denom)^(-1)
    return(q)
}
