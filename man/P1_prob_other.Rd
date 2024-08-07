% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/P1.R
\name{P1_prob_other}
\alias{P1_prob_other}
\title{Persistence probability via other approximations}
\usage{
P1_prob_other(
  R0,
  epsilon,
  burnout_prob_fun,
  k = 1,
  N = 10^6,
  subdivisions = 1000L,
  tiny = 0,
  ...
)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{epsilon}{mean infectious period (\eqn{1/(\gamma+\mu)}) as a
proportion of mean lifetime (\eqn{1/\mu}; \eqn{\varepsilon =
    0} corresponds to infinite lifetime, i.e., no mortality)}

\item{burnout_prob_fun}{function with which to calculate burnout
probability (either \code{\link{burnout_prob_vanH}} or
\code{\link{burnout_prob_MS}})}

\item{k}{initial number of infected individuals (\eqn{k})}

\item{N}{population size}

\item{subdivisions}{the maximum number of subintervals.}

\item{tiny}{avoid integration limits by this amount}

\item{...}{additional arguments pass to
\code{\link[stats]{integrate}}}
}
\description{
This computes \eqn{{\cal P}_1} using the burnout probability
approximations of van Herwaarden (1997) or Meerson and Sasorov
(2009).
}
\seealso{
\code{\link{burnout_prob_MS}},
\code{\link{burnout_prob_vanH}}, \code{\link{P1_prob}}
}
