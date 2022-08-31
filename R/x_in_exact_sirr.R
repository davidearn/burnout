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
