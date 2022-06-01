##' Plot persistence probability
##'
##' @inheritParams P1_prob
##' @param ... additional parameters passed to \code{\link{plot}}
##'
##' @export
##'
##' @examples
##' plot_P1(epsilon=0.05, N=10^6)
##' 
plot_P1 <- function(epsilon = 0.01, N = 10^6, npts = 1001,
                    P1_fun = P1_prob, add=FALSE, col="blue",
                    lwd=2, log="x", ... ) {
    xvals <- exp(seq(log(1.01), log(64), length=npts))
    P1vals <- P1_fun(xvals,epsilon=epsilon, N=N)
    if (any(Im(P1vals) != 0)) message("results are complex")
    print(summary(Im(P1vals)))
    cat("max(Im(P1)) = ", max(Im(P1vals)), "\n")
    P1vals <- Re(P1vals)
    if (!add) {
        plot(xvals, P1vals, las=1, col=col, type="l", lwd=lwd, log=log,
             xlab = "R0", ylab = "persistence probability", ...)
    } else {
        lines(xvals, P1vals, las=1, col="blue", type="l", lwd=lwd, ...)
    }
    legend("topleft",
           legend = c(sprintf("epsilon = %g", epsilon),
                      sprintf("N = %g", N)),
           bty="n")
    cat("max(Re(P1)) = ", max(P1vals), "\n")
}
