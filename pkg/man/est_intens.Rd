\name{est_intens}
\alias{est_intens}
\title{Estimates the intensity of the point pattern.}
\usage{
  est_intens(ppdata, covmatrix = NULL)
}
\arguments{
  \item{ppdata}{data of class ppp}

  \item{covmatrix}{(Optional) Covariance matrix of the
  kernel of a normal distribution.}
}
\value{
  A list of \item{ intensest }{ Estimated intensity (object
  of class "im", see \code{\link[spatstat]{density.ppp}}).
  } \item{ covmatrix }{ Covariance matrix. If
  \code{covmatrix = NULL} the matrix is estimated by
  \code{\link[ks]{Hscv}}. }
}
\description{
  Estimates the intensity of the point pattern by a kernel
  method (See \code{\link[spatstat]{density.ppp}}).
}
\examples{
data(craterA)
#change npixel = 50 with 1000 and you get a nicer picture
spatstat.options(npixel=50)
int <- est_intens(craterA)
plot(int$intensest, main = "pixel image of intensity")
plot(craterA$window, main = "contour plot of intensity"); contour(int$intensest, add =TRUE)
}
\seealso{
  \code{\link[spatstat]{density.ppp}},
  \code{\link[ks]{Hscv}}, \code{\link[spatstat]{eval.im}}
}

