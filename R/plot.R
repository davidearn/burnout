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

##' Plot comparison of two functions.
##'
##' \code{\link{plot}} method for \code{\link{compare_funs}} objects.
##' Typical use: exact vs approximate \eqn{q}.
##'
##' @inheritParams plot_P1
##' @param x \code{\link{compare_funs}} object, i.e., a
##'     \code{\link{data.frame}} with particular structure
##'
##' @seealso \code{\link{compare_funs}}, \code{\link{q_exact}},
##'     \code{\link{q_approx}}, \code{\link{x_in}}
##'
##' @import dplyr
##' @import ggplot2
##' 
##' @export
plot.compare_funs <- function(x, ...) {

    ## FIX: why am I getting this note from devtools::check():
    ##      plot.compare_q: no visible binding for global variable ‘R0’
    ##      I tried with and without pipe in case that somehow was the
    ##      but it makes no difference.

    ## scatter plot coloured by epsilon:
    ##scatter.plot <- (x %>% ggplot()
    scatter.plot <- (ggplot(data = x, ...)
        + geom_point(aes(x=f1, y=f2, colour=epsilon),
                     size = 0.5, alpha=0.5)
    )

    ##line.plot <- (x %>% ggplot()
    line.plot <- (ggplot(data = x, ...)
        + geom_point(aes(x=R0, y=f2, colour=epsilon))
        + geom_point(aes(x=R0, y=f1), colour = "red", size=0.25)
        + facet_wrap(~epsilon, scales="free_y")
        + scale_x_continuous(trans='log2')
    )

    ##relative.plot <- (x %>% ggplot()
    relative.plot <- (ggplot(data = x, ...)
        + geom_point(aes(x=R0, y=(f2-f1)/f1, colour=epsilon))
        + facet_wrap(~epsilon, scales="free_y")
        + scale_x_continuous(trans='log2')
    )

    print(scatter.plot)
    print(line.plot)
    print(relative.plot)
    return(invisible(list(scatter.plot,line.plot,relative.plot)))
}
