##' Plot persistence probability
##'
##' @inheritParams P1_prob
##' @param Rmin,Rmax minimum and maximum value of \eqn{{\cal R}_0}
##' @param show.N,show.epsilon if \code{TRUE}, display parameter value in legend
##' @param P1_fun function that computes \eqn{{\mathscr P}_1}
##' @param add if \code{TRUE}, add to existing plot
##' @param col,lwd,log standard \link{graphical parameters}
##' @param ... additional parameters passed to \code{\link{plot}}
##'
##' @importFrom graphics axis
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom latex2exp TeX
##'
##' @seealso P1_prob
##' 
##' @export
##'
##' @examples
##' op <- par(mfrow = c(2,2))
##' plot_P1(epsilon=0.01, N=10^4)
##' plot_P1(epsilon=0.01, N=10^5)
##' plot_P1(epsilon=0.01, N=10^6)
##' plot_P1(epsilon=0.01, N=10^7)
##' par(mfrow = op)
##' 
plot_P1 <- function(epsilon = 0.01, N = 10^6,
                    ##npts = 1001,
                    Rmin = 1.001,
                    Rmax = 64,
                    P1_fun = P1_prob, add=FALSE,
                    col="black", lwd=2, log="x",
                    show.N = TRUE,
                    show.epsilon = FALSE,
                    ... ) {

    ## naive R0 values:
    ##xvals <- exp(seq(log(Rmin), log(Rmax), length=npts))

    ## sprinkle R0 values where they are needed to get a smooth curve
    Rseq <- c(seq(Rmin,1.04,length=200),
              seq(1.04,1.1,length=100),
              seq(1.1,2,length=100),
              seq(2,Rmax,length=1000))
    
    P1vals <- P1_fun(Rseq, epsilon=epsilon, N=N)
    if (any(Im(P1vals) != 0)) {
        warning("results are complex")
        print(summary(Im(P1vals)))
        message("max(Im(P1)) = ", max(Im(P1vals)), "\n")
        P1vals <- Re(P1vals)
        cat("max(Re(P1)) = ", max(P1vals), "\n")
    }
    if (add) {
        lines(Rseq, P1vals, las=1, col=col, type="l", lwd=lwd,
              ...)
    } else {
        plot(Rseq, P1vals, las=1, col=col, type="l", lwd=lwd, log=log,
             bty="L", xaxs="i", xaxt="n",
             xlim=c(2^(-1/8),Rmax), ylim=c(10^(-4.5),1),
             xlab = expression(R[0]),
             ylab = "persistence probability",
             ...)
        xticks <- 2^(0:6)
        axis(side=1, at=xticks, labels=xticks)
        legend.text <- c()
        N.string <- sprintf("$N = 10^{%g}$", log10(N))
        epsilon.string <- sprintf("$\\varepsilon = %g$", epsilon)
        if (show.N)
            legend.text <- c(legend.text, latex2exp::TeX(N.string))
        if (show.epsilon)
            legend.text <- c(legend.text, latex2exp::TeX(epsilon.string))
        legend("topleft", legend = legend.text, bty="n")
    }
}
