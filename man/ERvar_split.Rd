% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ERvar_split.R
\name{ERvar_split}
\alias{ERvar_split}
\title{Extremely randomized split}
\usage{
ERvar_split(
  X,
  Y,
  ntry = 3,
  timeScale = 0.1,
  d_out = 0.1,
  FrechetSumOrMax = "max",
  ...
)
}
\arguments{
\item{X}{Inputs data}

\item{Y}{[list]: A list that contains the output, It must contain the
following elements (no choice): \code{type} defines the nature of the
output, can be "\code{curve}", "\code{sclalar}", "\code{factor}",
"\code{shape}", "\code{image}"; \code{Y} is the output variable; \code{id}
is the vector of the identifiers for each individuals, they should be the
same as the identifiers of the inputs.}

\item{ntry}{[numeric]: Only with \code{ERT=TRUE}, allows to manage with
randomness of the trees.}

\item{timeScale}{[numeric]: Allow to modify the time scale, increasing or
decreasing the cost of the horizontal shift. If timeScale is very big, then
the Frechet mean tends to the Euclidean distance. If timeScale is very
small, then it tends to the Dynamic Time Warping. Only used when there are
trajectories either in input or output.}

\item{d_out}{[numeric]: Time scale for the output curves (\code{d_out=0.1} by
default).}

\item{FrechetSumOrMax}{Frechet Mean and Frechet Distance can be define using
the 'sum' function or the 'max' function.}

\item{...}{: optional parameters to be passed to the low level function}
}
\description{
Extremely randomized split
}
\keyword{internal}
