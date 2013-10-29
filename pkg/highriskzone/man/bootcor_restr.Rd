\name{bootcor_restr}
\alias{bootcor_restr}
\title{Bootstrap correction to obtain desired failure probability}
\usage{
  bootcor_restr(ppdata, cutoff, numit = 100, tol = 0.001,
    nxprob = 0.1, hole = NULL, obsprobimage = NULL,
    intens = NULL, covmatrix = NULL, simulate = "intens",
    radiusClust = NULL, clustering = 5, verbose = TRUE)
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

  \item{hole}{(optional) an object of class \code{owin}
  representing a region inside the observation window of
  the \code{ppdata} where no observations were possible.}

  \item{obsprobsimage}{(optional) an object of class
  \code{im} giving the observation probabilities inside the
  observation window. Ranges of the coordinates must equal
  those of \code{ppdata}.  Only used if \code{obsprobs} is
  not given.}

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
\seealso{
  \code{\link[highriskzone]{det_hrz}},
  \code{\link[highriskzone]{eval_method}}

  The function offers the possibility to take into account
  so-called restriction areas. This is relevant in
  situations where the observed point pattern \code{ppdata}
  is incomplete. If it is known that no observations can be
  made in a certain area (for example because of water
  expanses), this can be accounted for by integrating a
  hole in the observation window. The shape and location of
  the hole is given by \code{hole}. Holes are part of the
  resulting high-risk zone. Another approach consists in
  weighting the observed events with their reciprocal
  observation probability when estimating the intensity. To
  do so, the observation probability can be specified by
  using \code{obsprobsimage} (image of the observation
  probability). Note that the observation probability may
  vary in space.

  For further information, see Mahling (2013), Appendix A
  (References).
}

