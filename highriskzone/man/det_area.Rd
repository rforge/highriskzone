\name{det_area}
\alias{det_area}
\title{Calculation of the area of the high-risk zone.}
\usage{
  det_area(win, intensmatrix, threshold)
}
\arguments{
  \item{win}{observation window}

  \item{intensmatrix}{matrix of the estimated intensity of
  the observed process (\code{as.matrix(intens)})}

  \item{threshold}{threshold c: The high-risk zone is the
  field in which the estimated intensity exceeds this
  value}
}
\value{
  A numerical value giving the area of the high-risk zone.
}
\description{
  This function is for the intensity-based method.
  Calculation of the area of the high-risk zone having the
  observation window, the intensity matrix and the
  threshold c. Used in function det_thresholdfromarea.
}
\seealso{
  \code{\link[spatstat]{owin}},
  \code{\link[spatstat]{area.owin}}
}

