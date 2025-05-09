% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak_prev.R
\name{peak_prev}
\alias{peak_prev}
\title{Peak prevalence}
\usage{
peak_prev(R0, epsilon = 0, eta.i = 0, eta.a = 0, xi = 1, yi = 0)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{epsilon}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean lifetime (\eqn{1/\mu}; \eqn{\varepsilon =
    0} corresponds to infinite lifetime, i.e., no mortality)}

\item{xi}{initial susceptible proportion \eqn{x_{\rm i}}}

\item{yi}{initial infective proportion \eqn{y_{\rm i}}}

\item{eta}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean duration of immunity (\eqn{1/\delta};
\eqn{\eta = 0} corresponds to permanent immunity}

\item{alpha}{antigenic evolution enhancement factor for effective
immunity waning; rate of immune decay is \eqn{\delta(1+\alpha
    Y)}.}
}
\value{
real number between 0 and 1
}
\description{
Approximation of peak prevalence in the SIR model with vital dynamics.
}
\details{
The peak prevalence is

\deqn{
  y_{\rm max} \approx 1-
    \frac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big)
    - \varepsilon
  \Big(1-\frac{1}{{\cal R}_0}\Big)
  \frac{(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) - 1}
  {(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) + {\cal R}_0} \,.
}

The expression is exact for \eqn{\varepsilon = 0}.
}
\examples{
peak_prev(R0=2)
peak_prev(R0=2, epsilon=c(0, 0.001, 0.01, 0.1, 0.9))
peak_prev(R0=c(2,4), epsilon=0.1)
peak_prev(R0=2, epsilon=0.001, eta=c(0.01, 0.1))
peak_prev(R0=c(2,4), epsilon=c(0,0.1)) # FIX: fails when both are vectors

}
\seealso{
\code{\link{peak_prev_nvd}}
}
