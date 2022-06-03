##' Susceptibles at boundary layer (approximation)
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
##'     \code{\link[expint]{expint_E1}}
##'
##' @inheritParams peak_prev
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
x_in <- function(R0, epsilon) {
    yeqm <- epsilon*(1-1/R0)
    ymax <- peak_prev(R0, epsilon)
    W0 <- emdbook::lambertW
    E1 <- expint::expint_E1
    xin <- -(1/R0)*W0(-R0*exp(R0*(yeqm-1)), b=0) +
        epsilon*exp(R0*yeqm)*(E1(R0*yeqm) - E1(R0*ymax))
    return(xin)
}
##
## from original code for Xb_approx in Xb.R:
##  Xb <- -(1/R)*lambertW(-R*exp(R*(Yb-1)), b=0) +
##               eps*exp(R*Yb)*(expint_E1(R*Yb) - expint_E1(R*Ymax))

##' Susceptibles at boundary layer (exact) [scalar version]
##'
##' @inheritParams x_in
##'
##' @importFrom sirr create_SIRmodel
##' @importFrom dplyr filter
##'
x_in_exact_scalar <- function(R0, epsilon, I0=1e-6) {
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
    mts <- compute_SIRts(m, saveRoots=TRUE, inits=c(S=1-I0), tmax=tmax)
    dd <- attr(mts, "rootPoints")
    Ihat.pts <- dd %>% filter(type == "Ieqm")
    xin <- Ihat.pts[2,"S"] # 2nd crossing
    ##Yb <- Ihat.pts[1,"I"]
    ##out <- c(Xb=Xb, Yb=Yb)
    return(xin)
}

##' Susceptibles at boundary layer (exact) [vector version]
##'
##' @inheritParams x_in
##'
##' @importFrom sirr create_SIRmodel
##' @importFrom dplyr filter
##'
##' @export
##' 
x_in_exact <- Vectorize(x_in_exact_scalar)
