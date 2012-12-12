\name{det_radius}
\alias{det_radius}
\title{Determination of the next neighbour distance which results in a high-risk zone with desired area}
\usage{
  det_radius(ppdata, distancemap, areahrz, win)
}
\arguments{
  \item{ppdata}{observed spacial point process of class
  ppp.}

  \item{distancemap}{distance map (object of class "im",
  see \code{\link[spatstat]{distmap}}): distance of every
  location in the observation window to the nearest other
  event}

  \item{areahrz}{given area of the high-risk zone}

  \item{win}{observation window of class owin}
}
\value{
  A list of \item{ cutoffdist }{ quantile of the next
  neighbour distance } \item{ thresh }{ distance }
}
\description{
  Used in function det_hrz
}
\seealso{
  \code{\link{det_nnarea}}, \code{\link[stats]{quantile}},
  \code{\link[stats]{uniroot}}
}

