% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/P1.R
\name{burnout_prob_hocb}
\alias{burnout_prob_hocb}
\title{Burnout probability using higher order corner/boundary
approximation of \eqn{x_{\rm in}}}
\usage{
burnout_prob_hocb(
  R0,
  epsilon,
  N = 10^6,
  xin = x_in_hocb(R0, epsilon),
  q_fun = q_approx
)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{epsilon}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean lifetime (\eqn{1/\mu}; \eqn{\varepsilon =
    0} corresponds to infinite lifetime, i.e., no mortality)}

\item{N}{population size}

\item{xin}{by default set via the approximation coded in \code{\link{x_in}}}

\item{q_fun}{function to compute Kendall's \eqn{q}}
}
\value{
real number between 0 and 1
}
\description{
This is just a convenient wrapper of \code{\link{burnout_prob}} to
use \code{\link{x_in_hocb}} rather than \code{\link{x_in}}.
}
\references{
\insertRef{Kendall1948b}{burnout}
}
\seealso{
\code{\link{burnout_prob}}
}
