\name{thin}
\alias{thin}
\title{Thinning of the observations (for evaluating the method)}
\usage{
  thin(full, nxprob)
}
\arguments{
  \item{full}{all observations of the point pattern}

  \item{nxprob}{probability of having unobserved events}
}
\value{
  A list of of observed and unobserved point patterns. Both
  of class ppp.
}
\description{
  The tinning is done by drawing independently from a
  Bernoulli distribution. Function is needed for functions
  eval_method, sim_clintens, sim_intens
}
\seealso{
  \code{\link[stats]{rbinom}}, \code{\link[spatstat]{ppp}}
  data(craterB) thdata <- thin(craterB, nxprob=0.1) thdata
  plot(thdata$observed); points(thdata$unobserved, col=4)
}

