% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fizzle.R
\name{fizzle_time}
\alias{fizzle_time}
\alias{t_delta}
\title{Fizzle time}
\usage{
fizzle_time(R0, k, delta)
}
\arguments{
\item{R0}{basic reproduction number (\eqn{{\cal R}_0})}

\item{k}{initial number of infected individuals (\eqn{k})}

\item{delta}{probability threshold}
}
\description{
Time \eqn{t_\delta} for which probability of fizzle after
\eqn{t_\delta} is \eqn{< \delta}
}
\details{
\deqn{
	t_{\delta} =  \frac{1}{{\cal R}_{0}-1} \ln\left(\frac{(1-\delta)^{-\frac{1}{k}}-\frac{1}{{\cal R}_{0}}}{(1-\delta)^{-\frac{1}{k}}-1}\right).
}
}
\examples{
fizzle_time(R0 = 2, k = 1, delta = 0.01)
plot_fizzle_time()

}
\seealso{
plot_fizzle_time
}
