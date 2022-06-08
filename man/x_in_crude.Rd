% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_in.R
\name{x_in_crude}
\alias{x_in_crude}
\title{Susceptibles at boundary layer (crude approximation)}
\usage{
x_in_crude(R0, epsilon, peakprev_fun = NULL, maxiter = 100)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{epsilon}{mean infectious period as a proportion of mean
lifetime (\eqn{\varepsilon = 0} corresponds to infinite
lifetime, i.e., no mortality)}

\item{peakprev_fun}{function of \eqn{{\cal R}_0} and
\eqn{\varepsilon} to use to compute peak prevalence}

\item{maxiter}{maximum numbers of iterations for convergence of
\code{\link[emdbook]{lambertW}}}
}
\description{
Susceptibles at boundary layer (crude approximation)
}
\details{
Argument \code{peakprev_fun} is ignored but is present so that the argument
list is the same as for \code{\link{x_in}}.
}