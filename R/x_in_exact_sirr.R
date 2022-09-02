
##' Susceptibles at boundary layer (exact) (scalar version)
##'
##' @export
##' @inheritParams x_in
##' @param I0 initial prevalence (proportion)
##'
x_in_exact_scalar <- function(R0, epsilon, I0=1e-6, tmax = NULL,
                               nt = 100,
                               return_traj = FALSE) {

    ## reasonable guess at an approximate max time for a specified R0
    ## (given default I0)
    if (is.null(tmax)) tmax <- 100/(R0-1)

    ## precompute useful factors, as well as hard-coding vars[i]
    ## lookup below rather than using with(as.list(...)), for speed
    herdImm <- if (R0 > 1) (1 - 1/R0) else 0

    ## epsilon == mu/(gamma + mu)
    ## -> epsilon*(gamma+mu) == mu
    ## -> gamma*epsilon = mu*(1-epsilon)
    ## -> mu = gamma*epsilon/(1-epsilon)
    Ihat <- epsilon * herdImm
    gamma <- 1
    mu <- gamma*epsilon/(1-epsilon)
    beta <- R0*(gamma+mu)

    sirgrad <- function(tau, vars, parms) {
        S <- vars[1]
        I <- exp(vars[2])
        dS <- mu * (1-S) - beta*S*I
        dlogI <- beta*S - (gamma + mu)
        return(list(c(dS=dS, dlogI=dlogI)))
    }

    eventfun <- function(t, x, p) { ..nroot <<- ..nroot + 1; return(x) }
    ## could also condition eventfun on (I = I^* AND deriv is negative)
    rfunc <- function(tau, vars, parms) {
        ## criterion for I == Ihat
        Ieqm <- exp(vars[2]) - Ihat
        ## eventfun is always called once initially (why?), so we want
        ## to make sure this criterion is non-zero until after the first root of Ieqm
        Ieqm2 <- if (..nroot <= 1) -1 else Ieqm
        return(c(Ieqm, Ieqm2))
    }

    ..nroot <- 0
    dd <- deSolve::lsodar(y = c(S=1-I0, logI = log(I0)),
                          ## ?? what should length.out be?
                          times = seq(0, tmax, length.out = nt),
                          func = sirgrad,
                          parms = list(R0 = R0, epsilon = epsilon),
                          rootfunc = rfunc,
                          events = list(func = eventfun, root = TRUE,
                                        ## terminate on criterion 2 (-> I==Ihat for second time)
                                        terminalroot = 2))

    if (return_traj) return(dd)
    tail(dd[,"S"], 1)
}

##'
##' Susceptibles at boundary layer (exact) (vector version)
##'
##' @rdname x_in_exact_scalar
##' @export
##' 
x_in_exact <- Vectorize(x_in_exact_scalar)


