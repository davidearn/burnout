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
##' \left(1+\dfrac{\varepsilon}{e^{\frac{{\cal R}_{0}}{\varepsilon} (1-x_{\rm in})}
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
q_exact <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    a <- (R0/epsilon)*(1-1/R0)
    x <- (R0/epsilon)*(1-xin)
    ##denom <- exp(x) * x^(-a) * lig(a,x)
    log.denom <- x - a*x + llig(a,x)
    denom <- exp(log.denom)
    q <- 1 / (1 + epsilon/denom)
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
##' \dfrac{1}{1+\sqrt{\frac{2\pi}{\varepsilon ({\cal R}_{0}-1)}} \;e^{\frac{{\cal R}_{0}}{\varepsilon}\big(\frac{1}{{\cal R}_{0}}-x_{\rm in}\big)}
##' \Big(\frac{1-\frac{1}{{\cal R}_{0}}}{1-x_{\rm in}}\Big)^{\frac{{\cal R}_{0}}{\varepsilon}\big(1-\frac{1}{{\cal R}_{0}}\big)}}\,.
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
q_approx <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    a <- (R0/epsilon)*(1-1/R0)
    b <- (R0/epsilon)*(1/R0-xin)
    ratio <- (1-1/R0) / (1-xin)
    log.messy <- (1/2)*(log(2*pi) - log(epsilon*(R0-1))) - b + a*log(ratio)
    denom <- 1 + exp(log.messy)
    q <- 1/denom
    return(q)
}

##' Compare exact and approximate calculations of Kendall's \eqn{q}
##'
##' @inheritParams peak_prev
##' @param n number of values of \code{R0} and \code{epsilon}
##' @param show.progress logical: if \code{TRUE} then spew each
##'     element of the data frame to \code{stdout} as it is computed
##'
##' @seealso \code{\link{q_exact}}, \code{\link{q_approx}}
##'
##' @export
##'
##' @return data frame with \eqn{\code{n}^2} rows
##'
##' @examples
##' (cq <- compare_q(n=3))
compare_q <- function(n=10, R0=seq(1.01,4,length=n),
                      epsilon=seq(0.0001,0.01,length=n),
                      show.progress=TRUE) {
    raw <- expand.grid(R0, epsilon)
    names(raw) <- c("R0", "epsilon")
    ## FIX: I obviously don't know how to use apply() ...
    ##      I should just use dplyr::mutate
    ## approx_fun <- function(x) q_approx(x['R0'], x['epsilon'])
    ## qapprox <- apply(raw, 1, approx_fun)
    ## print(qapprox) # FIX: debugging
    ## exact_fun <- function(x) q_exact(x['R0'], x['epsilon'])
    ## qexact <- apply(raw, 1, exact_fun)
    ## print(qexact) # FIX: debugging
    nn <- nrow(raw)
    qexact <- qapprox <- rep(NA,nn)
    if (show.progress) cat(sprintf("%s\t%s\t%s\t%s\t%s\n",
                                   "i", "R0", "epsilon", "qapprox", "qexact"))
    for (i in 1:nn) {
        r <- raw[i,"R0"]
        e <- raw[i,"epsilon"]
        if (show.progress) cat(sprintf("%d\t%g\t%g", i, r, e))
        qa <- q_approx(R0=r, epsilon=e)
        if (show.progress) cat(sprintf("\t%g", qa))
        qe <- q_exact(R0=r, epsilon=e)
        if (show.progress) cat(sprintf("\t%g\n", qe))
        qapprox[i] <- qa
        qexact[i] <- qe
    }
    dd <- cbind(raw, qexact, qapprox)
    return(dd)
}
