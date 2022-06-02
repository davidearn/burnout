
##' Duration of herd immunity
##'
##' The time \eqn{t_{\rm H}} (after the first epidemic) from when the
##' proportion susceptible drops below the herd immunity threshold
##' (defined via \eqn{X(t_{\rm H}) = 1/{\cal R}_0}) until the
##' proportion susceptible next rises above the threshold.
##'
##' \deqn{
##' 	t_{\rm H} =  \dfrac{1}{\varepsilon}
##'     \ln\Big( \dfrac{1-x_{\rm in}}{1 - \frac{1}{{\cal R}_0}} \Big)
##' }
##'
##' @inheritParams P1_prob
##' @inheritParams fizzle_prob
##' @param delta probability threshold
##'
##' @examples
##' op <- par(mfrow = c(1,2))
##' Rseq <- exp(seq(log(1.001),log(2),length=1001))
##' plot(Rseq, t_herd(Rseq,epsilon=0.01), type="l", lwd=2, bty="L",
##'      las=1, log="x", col="darkred", ylim=c(0,60),
##'      xlab = expression(R[0]), ylab = expression(t[H]))
##' lines(Rseq, t_herd(Rseq,epsilon=0.02), lwd=2, col="darkgreen")
##' lines(Rseq, t_herd(Rseq,epsilon=0.03), lwd=2, col="darkblue")
##' epsilonseq <- seq(0,0.03,length=1001)
##' plot(epsilonseq, t_herd(R0=1.5,epsilonseq), type="l", lwd=2, bty="L",
##'      las=1, col="darkred", ylim=c(0,60),
##'      xlab = expression(epsilon), ylab = expression(t[H]))
##' lines(epsilonseq, t_herd(R0=3,epsilonseq), lwd=2, col="darkgreen")
##' lines(epsilonseq, t_herd(R0=4.5,epsilonseq), lwd=2, col="darkblue")
##' par(mfrow = op)
##'
##' @export
##' 
t_herd <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    return((1/epsilon) * log((1-xin)/(1-1/R0)))
}
