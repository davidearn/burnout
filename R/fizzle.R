##' Fizzle probability
##' 
##' Probability that a disease will go extinct before causing a
##' substantial outbreak.
##'
##' The fizzle probability is
##' 
##' \deqn{1 - p_k =
##'   \begin{cases}
##'      1\,, & 0 \le {\cal R}_0 \le 1, \\
##'      \left(\dfrac{1}{{\cal R}_0}\right)^{k} \,, & 1\le{\cal R}_0.
##'   \end{cases}
##' }
##'
##' @param R0 basic reproduction number (\eqn{{\cal R}_0})
##' @param k initial number of infected individuals (\eqn{k})
##'
##' @return real number between 0 and 1
##' @export
##'
##' @seealso \code{\link{fizzle_prob}}
##'
##' @aliases fizz_prob
##' 
##' @examples
##' fizzle_prob(R0=c(1/2,2))
##' fizzle_prob(R0=2, k=2)
##' fizzle_prob(R0=c(2,4,8))
##' fizzle_prob(R0=2, k=c(1,2))  # FIX: fails with vector k
##' 
fizzle_prob <- function(R0, k=1) {
    ifelse( R0 < 1, 1, (1/R0)^k )
}

##' Not fizzle probability
##' 
##' Probability that a disease will \emph{not} go extinct before
##' causing a substantial outbreak.
##'
##' The probability of not fizzling is
##' 
##' \deqn{p_k =
##'   \begin{cases}
##'      0\,, & 0 \le {\cal R}_0 \le 1, \\
##'      1 - \left(\dfrac{1}{{\cal R}_0}\right)^{k} \,, & 1\le{\cal R}_0.
##'   \end{cases}
##' }
##'
##' @seealso \code{\link{fizzle_prob}}
##'
##' @param R0 basic reproduction number (\eqn{{\cal R}_0})
##' @param k initial number of infected individuals (\eqn{k})
##'
##' @return real number between 0 and 1
##' @export
##'
##' @aliases not_fizz_prob
##' 
##' @examples
##' not_fizzle_prob(R0=c(1/2,2))
##' not_fizzle_prob(R0=2, k=2)
##' not_fizzle_prob(R0=c(2,4,8))
##' not_fizzle_prob(R0=2, k=c(1,2))  # FIX: fails with vector k
##' 
not_fizzle_prob <- function(R0, k=1) {
    1 - fizzle_prob(R0, k)
}
