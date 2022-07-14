##' Susceptibles at boundary layer (boundary/boundary layer approximation)
##'
##' Fraction of hosts that are susceptible when the trajectory enters
##' the boundary layer at the end of the first epidemic.
##'
##' \deqn{
##' 	x_{\rm in} = -\dfrac{1}{{\cal R}_{0}} W_{0}\left(-{\cal R}_{0} e^{-{\cal R}_{0}\left(1-y^{\star}\right)}\right) + \varepsilon\, e^{{\cal R}_{0} y^{\star}}\big(E_{1}({\cal R}_{0} y^{\star}) - E_{1}({\cal R}_{0} y_{\rm max}) \big) \,.
##' }
##'
##' Here, \eqn{W_{0}} denotes the principal branch of the Lambert
##' \eqn{W} function
##' (\insertCite{Olver2010;textual}{burnout}, \insertCite{Corless1996;textual}{burnout},
##' [NIST 4.13](https://dlmf.nist.gov/4.13)),
##' \eqn{E_{1}(x) =
##' \int_{x}^{\infty} \frac{e^{-t}}{t}\,{\rm d}t} is the exponential
##' integral function (\insertCite{Olver2010;textual}{burnout}, [NIST
##' 6.2](https://dlmf.nist.gov/6.2)),
##' \eqn{y^\star=\varepsilon(1-\frac{1}{{\cal R}_0})}
##' is the equilibrium prevalence,
##' and \eqn{y_{\rm max}} is the \link[=peak_prev]{peak prevalence}.
##'
##' @seealso \code{\link{peak_prev}}, \code{\link[emdbook]{lambertW}},
##'     \code{\link[expint]{expint_E1}}, , \code{\link{x_in_cb}},
##'     \code{\link{x_in_exact}}
##'
##' @inheritParams peak_prev
##' @param peakprev_fun function of \eqn{{\cal R}_0} and
##'     \eqn{\varepsilon} to use to compute peak prevalence
##' @param maxiter maximum numbers of iterations for convergence of
##'     \code{\link[emdbook]{lambertW}}
##' @param ... additional arguments are ignored
##'
##' @return real number between 0 and 1
##' @importFrom emdbook lambertW
##' @importFrom expint expint_E1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Corless1996}{burnout}
##' 
##' \insertRef{Olver2010}{burnout}
##' 
##' @export
##'
##' @examples
##' x_in(R0=2, epsilon=0.01)
##' x_in(R0=2, epsilon=0)
##' x_in(R0=c(2,4), epsilon=0.01)
##' x_in(R0=2, epsilon=c(0.01, 0.1))
##' curve(x_in(x,epsilon=0.001), from=1.01, to=5, las=1, n=1001)
##' curve(x_in(x,epsilon=0.01), from=1.01, to=5, las=1, add=TRUE, col="magenta", n=1001)
##' curve(x_in(x,epsilon=0.1), from=1.01, to=5, las=1, add=TRUE, col="cyan", n=1001)
##'
x_in <- function(R0, epsilon, peakprev_fun = peak_prev, maxiter=100, ...) {
    yeqm <- epsilon*(1-1/R0) # equilibrium prevalence
    ymax <- peakprev_fun(R0, epsilon) # peak prevalence
    W0 <- function(x) {emdbook::lambertW(x, b=0, maxiter=maxiter, ...)}
    E1 <- expint::expint_E1
    xin <- -(1/R0)*W0(-R0*exp(R0*(yeqm-1))) +
        epsilon*exp(R0*yeqm)*(E1(R0*yeqm) - E1(R0*ymax))
    return(xin)
}
##
## from original code for Xb_approx in Xb.R:
##  Xb <- -(1/R)*lambertW(-R*exp(R*(Yb-1)), b=0) +
##               eps*exp(R*Yb)*(expint_E1(R*Yb) - expint_E1(R*Ymax))

##' Crude "Kermack-McKendrick" version of \eqn{x_{\rm in}}
##'
##' FIX: this is the same as \code{\link{x_in_crude}}, but coded
##' differently.
##'
##' This is just the case where \eqn{\varepsilon=0},
##' \deqn{x_{\rm in,KM}({\mathcal R}_0) = x_{\rm in}({\mathcal R}_0, 0)}
##'
##' @inheritParams x_in
##' @seealso \code{\link{x_in}}, \code{\link{x_in_crude}}
##'
##' @export
x_in_KM <- function(R0, maxiter=100, ...) {
    xin <- x_in(R0, epsilon=0, maxiter=maxiter, ...)
    return(xin)
}

##' Susceptibles at boundary layer (crude approximation)
##'
##' @details
##' Argument \code{peakprev_fun} is ignored but is present so that the argument
##' list is the same as for \code{\link{x_in}}.
##'
##' @inheritParams x_in
##' @importFrom emdbook lambertW
##' @export
##' 
x_in_crude <- function(R0, epsilon, peakprev_fun = NULL, maxiter=100, ...) {
    yeqm <- epsilon*(1-1/R0) # equilibrium prevalence
    W0 <- emdbook::lambertW
    xin <- -(1/R0)*W0(-R0*exp(R0*(yeqm-1)), b=0, maxiter=maxiter)
    return(xin)
}

##' Susceptibles at boundary layer (corner/boundary layer approximation)
##'
##' @details
##' \deqn{
##' 	x_{\rm in} = 1 + \left(1-\frac{1}{{\mathcal R}_{0}}\right) W_{-1}\left(-\dfrac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}} e^{-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}}\left(\dfrac{y_{\rm max}}{y^*}\right)^{\frac{\varepsilon}{{\mathcal R}_{0}}\frac{1}{1-\frac{1}{{\mathcal R}_{0}}}}\right)
##' }
##'
##' @inheritParams x_in
##' @importFrom emdbook lambertW
##' @seealso \code{\link{x_in}}, \code{\link{x_in_cb}}, \code{\link{x_in_exact}}
##' @export
##' 
x_in_cb <- function(R0, epsilon, peakprev_fun = peak_prev, maxiter=100, ...) {
    yeqm <- epsilon*(1-1/R0) # equilibrium prevalence
    ymax <- peakprev_fun(R0, epsilon) # peak prevalence
    W0 <- function(x) {emdbook::lambertW(x, b=0, maxiter=maxiter, ...)}
    Wm1 <- function(x) {emdbook::lambertW(x, b=-1, maxiter=maxiter, ...)}
    pc <- 1 - 1/R0 # p_crit
    xf <- -(1/R0)*W0(-R0*exp(-R0)) # standard final size
    Z <- 1 - xf
    xin <- 1 + pc * Wm1(-(Z/pc)*exp(-(Z/pc)) * (ymax/yeqm)^((epsilon/R0)/pc))
    return(xin)
}

##' \eqn{{\tilde Y}_1(x_{\rm f})}
##'
##' @details
##' \deqn{
##' \tilde{Y}_{1}(x{\rm f}) =  \int_{x_{\rm f}}^{1} \left(\frac{1}{{\mathcal R}_{0}t} -1\right) \frac{1-t}{{\mathcal R}_{0} t Y_{0}(t)} + \frac{1-x_{\rm f}}{{\mathcal R}_{0}x_{\rm f}}\frac{1}{t-x_{\rm f}} \, dt
##' }
##' 
##' @importFrom stats integrate
##'
##' @export
Ytilde_1 <- function(xf, R0, ...) {
    Y0 <- function(x) {1 - x - (1/R0)*log(x)}
    integrand <- function(t) {
        (1/(R0*t) - 1)*(1-t)/(R0*t*Y0(t)) + (1-xf)/(R0*xf)/(t-xf)
    }
    the.integral <-
        try(stats::integrate(f=integrand, lower=xf, upper=1,
                             ##subdivisions=subdivisions,
                             ...)$value)
    if ("try-error" %in% class(the.integral)) {
        warning("Ytilde_1: try error with xf = ", xf, "; returning NA")
        return(NA)
    }
    return(the.integral)
}

##' Susceptibles at boundary layer (higher order corner/boundary layer approximation)
##'
##' @details
##' \deqn{
##' 	x_{\rm in} = 1+\left(1-\frac{1}{{\mathcal R}_{0}}\right)
##'   W_{-1}\left(-\dfrac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}\left(\dfrac{1-x_{\rm f}}{y^*}\dfrac{\frac{1}{{\mathcal R}_{0}}-x_{\rm f}}{x_{\rm f}}\right)^{\frac{\varepsilon}{{\mathcal R}_{0}}\frac{1}{1-\frac{1}{{\mathcal R}_{0}}}}
##' e^{-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}-\varepsilon\frac{x_{\rm f}}{(1-x_{\rm f})\left(1-\frac{1}{{\mathcal R}_{0}}\right)}\tilde{Y}_{1}(x_{\rm f})}\right)
##' }
##'
##' @inheritParams x_in
##' @importFrom emdbook lambertW
##' @seealso \code{\link{x_in}}, \code{\link{x_in_cb}}, \code{\link{x_in_exact}}
##' @export
##' 
x_in_hocb <- function(R0, epsilon, peakprev_fun = peak_prev, maxiter=100, ...) {
    yeqm <- epsilon*(1-1/R0) # equilibrium prevalence
    ymax <- peakprev_fun(R0, epsilon) # peak prevalence
    W0 <- function(x) {emdbook::lambertW(x, b=0, maxiter=maxiter, ...)}
    Wm1 <- function(x) {emdbook::lambertW(x, b=-1, maxiter=maxiter, ...)}
    pc <- 1 - 1/R0 # p_crit
    xf <- -(1/R0)*W0(-R0*exp(-R0)) # standard final size
    Z <- 1 - xf
    xin <- 1 + pc * Wm1(-(Z/pc)*((Z/yeqm)*((1/R0-xf)/xf))^((epsilon/R0)/pc) * exp(-(Z/pc) - epsilon*(xf/(Z*pc))*Ytilde_1(xf,R0)))
    return(xin)
}

##' Susceptibles at boundary layer (exact) (scalar version)
##'
##' @inheritParams x_in
##' @param I0 initial prevalence (proportion)
##'
##' @importFrom sirr create_SIRmodel
##' @importFrom sirr compute_SIRts
##' @importFrom dplyr filter
##
## The following avoids check() from complaining "no visible binding
## for global variable ‘type’".
## https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887/2
## https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
##' @importFrom rlang .data
##'
x_in_exact_scalar <- function(R0, epsilon, I0=1e-6, ...) {
    m <- sirr::create_SIRmodel(R0=R0, eps=epsilon)
    ## FIX: tmax=100 works for a large range of R0, but is still too
    ## short as R0 --> 1.  I need to have a stopcrit that ends when
    ## the 2nd crossing occurs.
    tmax <- 100
    if (R0 < 1.5) tmax <- 200
    if (R0 < 1.1) tmax <- 300
    if (R0 < 1.05) tmax <- 1000
    if (R0 < 1.01) tmax <- 2000
    if (R0 < 1.005) tmax <- 20000
    mts <- sirr::compute_SIRts(m, saveRoots=TRUE, inits=c(S=1-I0), tmax=tmax)
    dd <- attr(mts, "rootPoints")
    Ihat.pts <- dd %>% filter(.data$type == "Ieqm")
    xin <- Ihat.pts[2,"S"] # 2nd crossing
    ##Yb <- Ihat.pts[1,"I"]
    ##out <- c(Xb=Xb, Yb=Yb)
    return(xin)
}

##' Susceptibles at boundary layer (exact) (vector version)
##'
##' @inheritParams x_in_exact_scalar
##'
##' @importFrom sirr create_SIRmodel
##' @importFrom dplyr filter
##'
##' @export
##' 
x_in_exact <- Vectorize(x_in_exact_scalar)


##' Plot \eqn{x_{\rm in}}
##' 
##' @inheritParams x_in
##' @inheritParams graphics::par
##' @inheritParams compare_funs
##'
##' @param xvar string indicating whether \code{R0} or \code{epsilon}
##'     is the desired independent variable
##' @param col,log,lwd,xlim,ylim,... see \code{\link{graphical parameters}}
##' @param xin_fun function of \code{R0} and \code{epsilon} that
##'     returns susceptible proportion at entry to boundary layer
##'
##' @seealso \code{\link{x_in}}, \code{\link{x_in_exact}}
##'
##' @importFrom graphics title
##'
##' @export
##'
##' @examples
##' op <- par(mfrow = c(1,2))
##' plot_x_in(epsilon = 0.01)
##' plot_x_in(epsilon = 0.1)
##' par(mfrow = op)
##' op <- par(mfrow = c(1,3))
##' plot_x_in(R0 = 2, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' plot_x_in(R0 = 5, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' plot_x_in(R0 = 20, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' par(mfrow = op)
##'
plot_x_in <- function(Rmin=1.0001, Rmax=20,
                      ##R0 = exp(seq(log(Rmin),log(Rmax),length=1001)),
                      R0 = seq(Rmin,Rmax,length=1001),
                      epsilon = c(0.01, 0.02, 0.03),
                      xvar = "R0"
                    , xin_fun = x_in
                    , col = 1:length(epsilon)
                      ##col = c("darkred", "darkgreen", "darkblue")
                    , lwd=2
                    , log="x"
                    , xlim = if (xvar == "R0") c(1,Rmax) else c(0,1)
                    , ylim = if (xvar == "R0") c(0,1) else c(0,max(1/R0))
                    , ...
                      ) {
    if (xvar == "R0") {
        cat("plot_x_in: Rmin = ", Rmin, ", Rmax = ", Rmax, "\n")
        xx <- R0
    } else {
        cat("plot_x_in: plotting as a function of epsilon...\n")
        xx <- epsilon
        if (length(R0) != 1) stop("only one R0 value allow if xvar is epsilon")
    }
    
    ## show naive approx (eqm susceptible proportion) as dashed line first:
    plot(xx,
         if (xvar == "R0") 1/R0 else rep(1/R0,length(xx)),
         type="l", log=log, lwd=lwd/2, bty="L", lty="dashed", las=1,
         xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
         xlab = if (xvar == "R0") expression(R[0]) else expression(epsilon),
         ylab = expression(x[i][n]), ...)
    title(main = "Susceptible proportion at\nentry to boundary layer")

    if (xvar == "R0") {
        ## curves with peak prevalence estimate from approximation of integral:
        for (iepsilon in 1:length(epsilon)) {
            ## estimate obtained by taking x_in to be x_f from final size formula:
            lines(R0, x_in_crude(R0,epsilon=epsilon[iepsilon]), lwd=lwd/2, col="pink")
            ## better estimate
            lines(R0, xin_fun(R0,epsilon=epsilon[iepsilon]), lwd=lwd, col=col[iepsilon])
        }
        ## plot with dotted curves when using the crude peak prevalence
        ## ignoring vital dynamics:
        for (iepsilon in 1:length(epsilon)) {
            lines(R0, xin_fun(R0,epsilon=epsilon[iepsilon], peakprev_fun=peak_prev_nvd),
                  lwd=lwd,
                  col=col[iepsilon],
                  lty="dotted")
        }
        legend("topright", bty="n", title=expression(epsilon),
               legend = epsilon, col = col, lwd=lwd)
    } else {
        ## estimate obtained by taking x_in to be x_f from final size formula:
        lines(xx, x_in_crude(R0,xx), lwd=lwd/2, col="pink")
        ## better estimate
        lines(xx, xin_fun(R0,xx), lwd=lwd, col=col)
        lines(xx, xin_fun(R0,xx, peakprev_fun=peak_prev_nvd),
              lwd=lwd, col=col, lty="dotted")
        title(sub = latex2exp::TeX(sprintf("$R_0 = %g$", R0)))
    }
}

##' Plot all approximations to \eqn{x_{\rm in}}
##' 
##' @inheritParams plot_x_in
##' @inheritParams graphics::par
##' @inheritParams compare_funs
##'
##' @seealso \code{\link{plot_x_in}}, \code{\link{x_in}},
##'     \code{\link{x_in_exact}}
##'
##' @importFrom graphics title
##'
##' @export
##'
plot_all_x_in <- function(Rmin=1.0001, Rmax=20,
                      ##R0 = exp(seq(log(Rmin),log(Rmax),length=1001)),
                      R0 = seq(Rmin,Rmax,length=1001),
                      epsilon = c(0.01, 0.02, 0.03),
                      xvar = "R0"
                    , col = 1:length(epsilon)
                      ##col = c("darkred", "darkgreen", "darkblue")
                    , lwd=2
                    , log="x"
                    , xlim = if (xvar == "R0") c(1,Rmax) else c(0,1)
                    , ylim = if (xvar == "R0") c(0,1) else c(0,max(1/R0))
                    , ...
                      ) {
    omf <- par(mfrow = c(2,2))
    on.exit(par(mfrow = omf))

    ##plot_x_in(xin_fun = x_in_exact) # too slow
    plot_x_in(xin_fun = x_in_crude)
    legend("top", bty="n", legend="x_in_crude")
    plot_x_in(xin_fun = x_in)
    legend("top", bty="n", legend="x_in")
    plot_x_in(xin_fun = x_in_cb)
    legend("top", bty="n", legend="x_in_cb")
    plot_x_in(xin_fun = x_in_hocb)
    legend("top", bty="n", legend="x_in_hocb")
}
