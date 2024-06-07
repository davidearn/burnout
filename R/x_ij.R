Euler_gamma <- -digamma(1)
Ein <- function( z ) {E1(z) + log(z) + Euler_gamma} 

Winv <- function( x ) {x * exp(x)}

## 

##' Recursive function to compute \eqn{x_{{\rm i},j}}
##'
##' This uses the recurrence relation from Parsons and Earn (2023).
##' FIX: cite paper.
##'
##' @seealso \code{\link{P1_prob}}
##' 
##' @inheritParams P1_prob
##' @param j epidemic wave index
##'
##' @export
##'
x_i <- function ( R0, j ) {
    if (j<1) stop("x_i: j = ", j)
    if (j==1) return( 1 ) # DFE
    ## now we know j>1:
    xstar <- 1/R0
    xijm1 <- x_i( R0, j-1 )
    xfjm1 <- x_f( R0, j-1 )
    return( 1 + (1-xstar) * W0(Winv(-(1-xfjm1)/(1-xstar))) )
}

x_f <- function( R0, j ) {
    xstar <- 1/R0
    xij <- x_i( R0, j )
    return( -xstar * W0( Winv(-xij/xstar) ) )
}
