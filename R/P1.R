##' Persistence probability
##'
##' Probability of surviving the trough after the first epidemic.
##'
##' \deqn{
##'     {{\mathscr P}_1}({\cal R}_0,\varepsilon,k,n) =
##'       p_k\big( 1 - q(x_{\rm in})^{ny^\star} \big)
##' }
##' where \eqn{k} is the initial number of infected individuals and \eqn{n}
##' is the population size.
##'
##' @seealso \code{\link{fizzle_prob}}, \code{\link{x_in}},
##'     \code{\link{llig}}, \code{\link{q_approx}}
##'
##' @inheritParams q_approx
##' @param n population size
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
##' P1_prob(R0=2, epsilon=0.01)
##' P1_prob(R0=2, epsilon=0)
##' P1_prob(R0=c(2,4), epsilon=0.01)
##' P1_prob(R0=2, epsilon=c(0.01, 0.1))
##' curve(P1_prob(x,epsilon=0.001), from=1.01, to=2, las=1, n=1001, ylim=c(0,1))
##' curve(P1_prob(x,epsilon=0.01), from=1.01, to=2, las=1, add=TRUE, col="magenta", n=1001)
##' curve(P1_prob(x,epsilon=0.1), from=1.01, to=2, las=1, add=TRUE, col="cyan", n=1001)
##' xvals <- exp(seq(log(1.01), log(64), length=1001))
##' P1vals <- P1_prob(xvals,epsilon=0.01)
##' ## results are complex:
##' max(Im(P1vals))
##' P1vals <- Re(P1vals)
##' plot(xvals, P1vals, las=1, col="blue", type="l", lwd=2, log="x")
##' cat("max(P1) = ", max(P1vals), "\n")
##'
P1_prob <- function(R0, epsilon, k=1, n=10^6, xin = x_in(R0,epsilon),
                    q_fun = q_approx) {
    pk <- 1 - fizzle_prob(R0, k)
    ystar <- epsilon * (1 - 1/R0)
    q <- q_fun(R0, epsilon, xin)
    P1 <- pk * (1 - q^(n*ystar))
    return(P1)
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