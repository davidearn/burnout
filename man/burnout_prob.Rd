% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/P1.R
\name{burnout_prob}
\alias{burnout_prob}
\title{Burnout probability based on Kendall's \eqn{q}}
\usage{
burnout_prob(
  R0,
  epsilon,
  eta = 0,
  alpha = 0,
  N = 10^6,
  xin = x_in(R0 = R0, epsilon = epsilon, eta = eta, alpha = alpha),
  q_fun = q_approx,
  ystar = eqm_prev(R0 = R0, epsilon = epsilon, eta = eta, alpha = alpha)
)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{epsilon}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean lifetime (\eqn{1/\mu}; \eqn{\varepsilon =
    0} corresponds to infinite lifetime, i.e., no mortality)}

\item{eta}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean duration of immunity (\eqn{1/\delta};
\eqn{\eta = 0} corresponds to permanent immunity}

\item{alpha}{antigenic evolution enhancement factor for effective
immunity waning; rate of immune decay is \eqn{\delta(1+\alpha
    Y)}.}

\item{N}{population size}

\item{xin}{by default set via the approximation coded in \code{\link{x_in}}}

\item{q_fun}{function to compute Kendall's \eqn{q}}

\item{ystar}{prevalence at boundary layer (equilibrium prevalence
by default)}
}
\value{
real number between 0 and 1
}
\description{
\deqn{
    q^{Ny^\star}
}
}
\references{
\insertRef{Kendall1948b}{burnout}
}
\seealso{
\code{\link{fizzle_prob}}, \code{\link{x_in}},
\code{\link{llig}}, \code{\link{q_approx}}, \code{\link{q_exact}}
}
