
##' Duration of herd immunity
##'
##' The time \eqn{t_{\rm H}} (after the first epidemic) from when the
##' proportion susceptible drops below the herd immunity threshold
##' (defined via \eqn{X(t_{\rm H}) = 1/{\cal R}_0}) until the
##' proportion susceptible next rises above the threshold.
##'
##' \deqn{
##' 	t_{\rm H} =  \frac{1}{\varepsilon}
##'     \ln\Big( \frac{1-x_{\rm in}}{1 - \frac{1}{{\cal R}_0}} \Big)
##' }
##'
##' @inheritParams P1_prob
##' @inheritParams fizzle_prob
##'
##' @examples
##' op <- par(mfrow = c(1,2))
##' plot_t_herd(col = c("darkred", "darkgreen", "darkblue"))
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
    ##return((1/epsilon) * log((1-xin)/(1-1/R0)))
    return((1/epsilon) * (log(1-xin) - log(1-1/R0)))
}

##' Plot herd immunity duration
##'
##' @inheritParams t_herd
##' @inheritParams base::plot
##' @param col,lwd,ylim,... see \code{\link{graphical parameters}}
##'
##' @seealso \code{\link{t_herd}}
##'
##' @importFrom graphics title
##'
##' @export
##'
##' @examples
##' plot_t_herd(delta = 0.01)
##'
plot_t_herd <- function(R0 = exp(seq(log(1.001),log(2),length=1001)),
                        epsilon = c(0.01, 0.02, 0.03)
                      , col = 1:length(epsilon)
                        ##col = c("darkred", "darkgreen", "darkblue")
                      , lwd=2
                        ##, log="x"
                      , ylim=c(0,60)
                      , ...
                        ) {
    plot(R0, t_herd(R0,epsilon=epsilon[1]), type="l", lwd=lwd, bty="L",
         las=1, col=col[1], ylim=ylim,
         xlab = expression(R[0]), ylab = expression(t[H]), ...)
    title(main = latex2exp::TeX("Duration of herd immunity $t_{H}$"))
    for (iepsilon in 2:length(epsilon)) {
        lines(R0, t_herd(R0,epsilon=epsilon[iepsilon]), lwd=2, col=col[iepsilon])
    }
    legend("topright", bty="n", title=expression(epsilon),
           legend = epsilon, col = col, lwd=lwd)
}
