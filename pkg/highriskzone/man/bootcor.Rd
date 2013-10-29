\name{bootcor}
\alias{bootcor}
\title{Bootstrap correction to obtain desired failure probability}
\usage{
  bootcor(ppdata, cutoff, numit = 100, tol = 0.001,
    nxprob = 0.1, intens = NULL, covmatrix = NULL,
    simulate = "intens", radiusClust = NULL,
    clustering = 5, verbose = TRUE)
}
\arguments{
  \item{ppdata}{Observed spatial point process of class
  ppp.}

  \item{cutoff}{Desired failure probability alpha, which is
  the probability of having unobserved events outside the
  high-risk zone.}

  \item{numit}{Number of iteration to perform (per tested
  value for cutoff)}

  \item{tol}{Tolerance: acceptable difference between the
  desired failure probability and the fraction of high-risk
  zones not covering all events}

  \item{nxprob}{Probability of having unobserved events.
  Default value is 0.1.}

  \item{intens}{(optional) estimated intensity of the
  observed process (object of class "im", see
  \code{\link[spatstat]{density.ppp}}), only needed for
  type="intens". If not given, it will be estimated.}

  \item{covmatrix}{(optional) Covariance matrix of the
  kernel of a normal distribution, only meaningful for
  \code{type="intens"} if no intensity is given. If not
  given, it will be estimated.}

  \item{simulate}{The type of simulation, can be one of
  \code{"thinning", "intens"} or \code{"clintens"}}

  \item{radiusClust}{(Optional) radius of the circles
  around the parent points in which the cluster points are
  located. Only used for \code{simulate = "clintens"}.}

  \item{clustering}{a value >= 1 which describes the amount
  of clustering; the adjusted estimated intensity of the
  observed pattern is divided by this value; it also is the
  parameter of the Poisson distribution for the number of
  points per cluster. Only used for \code{simulate =
  "clintens"}.}

  \item{verbose}{logical. Should information on tested
  values/progress be printed?}
}
\value{
  An object of class bootcorr, which consists of a list of
  the final value for alpha (alphastar) and a data.frame
  course containing information on the simulation course,
  e.g. the tested values.
}
\description{
  Description
}
\examples{
data(craterB)
set.seed(4321)

bc <- bootcor(ppdata=craterB, cutoff=0.2, numit=100, tol=0.02, nxprob=0.1)
bc
summary(bc)
plot(bc)
}
\seealso{
  \code{\link[highriskzone]{det_hrz}},
  \code{\link[highriskzone]{eval_method}}
}

