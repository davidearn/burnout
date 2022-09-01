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

## epsilon: mean generation interval as fraction of mean lifetime
## R0
## -> mean lifetime = L = 1/mu (= 1, nondimensional?)
## -> mean generation interval = epsilon*L = 1/gamma
## -> beta/(gamma+mu) = R0 -> beta = R0*(gamma+mu) = R0*(1/epsilon +1)
sirgrad <- function(tau, vars, parms) {
    with(as.list(c(parms, vars)), {
        ## more efficient to precompute these ...
        mu <- 1
        gamma <- 1/epsilon
        beta <- R0*(gamma+mu)
        I <- exp(logI)
        dS <- mu * (1-S) - beta*S*I
        dlogI <- beta*S - (gamma + mu)
        return(list(c(dS=dS, dlogI=dlogI)))
  })
}

## 1-1/R0 = I*(1/epsilon + 1)
## I = (1-1/R0)/(1/epsilon + 1) = (1-1/R0)*epsilon/(1+epsilon)
rfunc <- function(tau, vars, parms) {
    with(as.list(c(parms, vars)), {
        herdImm <- if (R0 > 1) (1 - 1/R0) else 0
        Ihat <- epsilon/(1+epsilon) * herdImm
        return(c(Ieqm = exp(logI) - Ihat))
    })
}

if (FALSE) {

    pfun <- function(d, parms = p0) {
        par(las = 1, bty = "l")
        plot(logI ~ time, data = as.data.frame(d), type ="l", log = "x")
        Ihat <- with(as.list(parms), epsilon/(1+epsilon)*(1-1/R0))
        abline(h = log(Ihat), col = 2)
    }


    ## shouldn't eventfun only be called at a root?
    eventfun <- function(t, x, p) { nroot <<- nroot + 1; return(x) }  ## no-op/unchanged
    ## 1-1/R0 = I*(1/epsilon + 1)
    ## I = (1-1/R0)/(1/epsilon + 1) = (1-1/R0)*epsilon/(1+epsilon)
    rfunc <- function(tau, vars, parms) {
        with(as.list(c(parms, vars)), {
            herdImm <- if (R0 > 1) (1 - 1/R0) else 0
            Ihat <- epsilon/(1+epsilon) * herdImm
            Ieqm <- exp(logI) - Ihat
            Ieqm2 <- if (nroot<2) -1 else Ieqm
            ## cat(Ieqm, nroot, Ieqm2, "\n")
            return(c(Ieqm, Ieqm2))
        })
    }

    I0 <- 1e-6
    p0 <- c(R0 = 5, epsilon = 0.1)
    d1 <- deSolve::ode(y = c(S=1-I0, logI = log(I0)),
                       times = seq(0, 100, length.out = 5001),
                       func = sirgrad,
                       parms = p0)

    pfun(d1)
    nroot <- 0

    d2 <- deSolve::lsodar(y = c(S=1-I0, logI = log(I0)),
                          times = seq(0, 10, length.out = 501),
                          func = sirgrad,
                          parms = list(R0 = 5, epsilon = 0.1),
                          rootfunc = rfunc,
                          events = list(func = eventfun, root = TRUE, terminalroot = 2))
    pracma::toc()
    
    pfun(d2)
    tail(d2[,"time"], 1)  ## time of second root

}
