##' Persistence probability
##'
##' Probability of surviving the trough after the first epidemic.
##'
##' If \eqn{p_k} is the probability of \emph{not}
##' \code{\link[=fizzle_prob]{fizzling}} if \eqn{k} individuals are
##' initially infected, then the probability of burning out
##' conditional on not fizzling is
##'
##' \deqn{ p_k q^{Ny^\star} \,,}
##'
##' where \eqn{q} is Kendall's \eqn{q}, the expected total population
##' size is \eqn{N}, and the equilibrium prevalence is
##' 
##' \deqn{y^\star = \varepsilon\Big(1 - \dfrac{1}{{\cal R}_0}\Big) \,.}
##' 
##' Hence, the probability of \emph{either} fizzling before a first
##' major outbreak \emph{or} burning out after a first major outbreak
##' is
##'
##' \deqn{
##'      (1-p_k) + p_k q^{Ny^\star} =
##'      1 - p_k\big( 1 - q^{Ny^\star} \big) \,.
##' }
##' 
##' Consequently, the probability of \emph{persisting} beyond the
##' first epidemic is
##' 
##' \deqn{
##'     {{\mathscr P}_1}({\cal R}_0,\varepsilon,k,N) =
##'       p_k\big( 1 - q^{Ny^\star} \big) \,.
##' }
##'
##' @seealso \code{\link{fizzle_prob}}, \code{\link{x_in}},
##'     \code{\link{llig}}, \code{\link{q_approx}}, \code{\link{q_exact}}
##'
##' @inheritParams q_approx
##' @inheritParams fizzle_prob
##' @param N population size
##' @param q_fun function to compute Kendall's \eqn{q}
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
P1_prob <- function(R0, epsilon, k=1, N=10^6, xin = x_in(R0,epsilon),
                    q_fun = q_approx) {
    fizz <- fizzle_prob(R0, k)
    ## pk = probability of not fizzling:
    notfizz <- 1 - fizz
    ystar <- epsilon * (1 - 1/R0)
    ## Kendall's q:
    q <- q_fun(R0, epsilon, xin)
    ## probability of not fizzling and then burning out:
    notfizz.and.burn <- notfizz * q^(N*ystar)
    ## probability of either fizzling or burning out:
    fizz.or.burn <- fizz + notfizz.and.burn
    ## persist after neither fizzling nor burning out:
    P1 <- 1 - fizz.or.burn # P1 <- pk * (1 - q^(N*ystar))
    return(P1)
}
