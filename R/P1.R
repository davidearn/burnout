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
