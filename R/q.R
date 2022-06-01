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
    log.denom <- x - a*log(x) + llig(a,x)
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
##' \left(1 +
##' \dfrac{1}{1+\sqrt{\frac{2\pi}{\varepsilon ({\cal R}_{0}-1)}} \;e^{\frac{{\cal R}_{0}}{\varepsilon}\big(\frac{1}{{\cal R}_{0}}-x_{\rm in}\big)}
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
##q_approx_ORIG <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
q_approx <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    a <- (R0/epsilon)*(1 - 1/R0)
    b <- (R0/epsilon)*(1/R0 - xin)
    log.fac1 <- (1/2)*(log(2*pi) - log(epsilon*(R0-1)))
    log.fac2 <- a*(log(1-1/R0) - log(1-xin))
    log.messy <- log.fac1 + log.fac2 + b
    denom <- 1 + exp(log.messy)
    #q <- 1/denom
    q <- (1 + 1/denom)^(-1)
    return(q)
}

##' Compare exact vs approximate calculations of Kendall's \eqn{q}
##'
##' @details The default \code{R0} value of \code{NULL} yields a
##'     sensible sprinkling of \code{nR0} \eqn{{\cal R}_0} values.
##'
##' @inheritParams peak_prev
##' @param Rmin,Rmax minimum and maximum values of \eqn{{\cal R}_0}
##' @param nR0,nepsilon number of values of \code{R0} and \code{epsilon}
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
##' (cq <- compare_q(nR0=3))
##' 
compare_q <- function(R0 = NULL, Rmin = 1.001, Rmax = 64, nR0 = 101,
                      epsilon = 10^(-(4:1)), nepsilon = length(epsilon),
                      show.progress=FALSE) {

    if (is.null(R0)) {
        ## sprinkle R0 values where they are needed to get a smooth curve
        n1 <- round(0.2*nR0)
        n2 <- round(0.1*nR0)
        R0 <- c(seq(Rmin,1.04,length=n1),
                seq(1.04,1.1,length=n2),
                seq(1.1,2,length=n2),
                exp(seq(log(2),log(8),length=n2)),
                seq(8,Rmax,length=(nR0-3*n2-n1)))
    } else {
        nR0 <- length(R0)
    }

    raw <- expand.grid(R0, epsilon)
    names(raw) <- c("R0", "epsilon")
    ## FIX: The following would be much cleaner if dplyr::mutate were
    ##      applied to the raw data frame, but the loop allows me to
    ##      print the table as far as it can be be computed, which is
    ##      helpful.
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
    ## FIX: not sure why this is failing...
    ## nIme <- sum(Im(qexact) != 0)
    ## nIma <- sum(Im(qapprox) != 0)
    ## if (nIme != 0) warning(nIme, " complex values of qexact")
    ## if (nIma != 0) warning(nIma, " complex values of qapprox")
    dd <- cbind(raw, Re(qexact), Re(qapprox))
    class(dd) <- c("compare_q", "data.frame")
    ## attr(dd,"nIme") <- nIme
    ## attr(dd,"nIma") <- nIma
    return(dd)
}
