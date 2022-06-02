
##' Duration of herd immunity
##'
##' The time \eqn{t_{\rm H}} (after the first epidemic) from when the
##' proportion susceptible drops below the herd immunity threshold
##' (defined via \eqn{X(t_{\rm H}) = 1/{\cal R}_0}) until the
##' proportion susceptible next rises above the threshold.
##'
##' \deqn{
##' 	t_{\rm H} =  \dfrac{1}{\varepsilon}
##'     \ln\Big( \dfrac{1-x_{\rm in}}{1 - \frac{1}{{\cal R}_0}} \Big)
##' }
##'
##' @inheritParams P1_prob
##' @inheritParams fizzle_prob
##' @param delta probability threshold
##'
##' @export
##' 
t_herd <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    return((1/epsilon) * log((1-xin)/(1-1/R0)))
}
