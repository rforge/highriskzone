#' Bootstrap correction to obtain desired failure probability
#'
#' Description
#'
#' @param ppdata Observed spatial point process of class ppp.
#' @param cutoff Desired failure probability alpha, which is the probability of having
#'                 unobserved events outside the high-risk zone.
#' @param numit Number of iteration to perform (per tested value for cutoff)
#' @param tol Tolerance: acceptable difference between the desired failure probability and the fraction of
#'             high-risk zones not covering all events
#' @param nxprob Probability of having unobserved events.
#'                Default value is 0.1.
#' @param hole  (optional) an object of class \code{owin} representing a region inside the observation window of
#'               the \code{ppdata} where no observations were possible.
#' @param obsprobsimage  (optional) an object of class \code{im} giving the observation probabilities inside the
#'                        observation window. Ranges of the coordinates must equal those of \code{ppdata}.
#'                        Only used if \code{obsprobs} is not given.
#' @param intens (optional) estimated intensity of the observed process (object of class "im",
#'                see \code{\link[spatstat]{density.ppp}}), only needed for type="intens". If not given,
#'                it will be estimated.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only meaningful for
#'                    \code{type="intens"} if no intensity is given. If not given, it will be estimated.
#' @param simulate The type of simulation, can be one of \code{"thinning", "intens"} or \code{"clintens"}
#' @param radiusClust (Optional) radius of the circles around the parent points in which the cluster
#'                    points are located. Only used for \code{simulate = "clintens"}.
#' @param clustering a value >= 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it also is the parameter of the Poisson distribution
#'          for the number of points per cluster. Only used for \code{simulate = "clintens"}.
#' @param verbose logical. Should information on tested values/progress be printed?
#' @expor
#' @return An object of class bootcorr, which consists of a list of the final value for alpha (alphastar)
#'         and a data.frame course containing information on the simulation course, e.g. the tested values.
#' @seealso \code{\link[highriskzone]{det_hrz}}, \code{\link[highriskzone]{eval_method}}
#'
#'
#' The function offers the possibility to take into account so-called restriction areas. This is relevant in
#' situations where the observed point pattern \code{ppdata} is incomplete. If it is known that no observations
#' can be made in a certain area (for example because of water expanses),
#' this can be accounted for by integrating a hole in the observation window.
#' The shape and location of the hole is given by \code{hole}. Holes are
#' part of the resulting high-risk zone.
#' Another approach consists in weighting the observed events with their reciprocal observation probability when
#' estimating the intensity. To do so, the observation probability can be specified by using
#' \code{obsprobsimage} (image of the observation probability). Note that the
#' observation probability may vary in space.
#'
#' For further information, see Mahling (2013), Appendix A (References).
#'


bootcor_restr <- function(ppdata, cutoff, numit = 100, tol=0.001, nxprob = 0.1,
                        hole=NULL, obsprobimage=NULL,
                        intens = NULL, covmatrix = NULL, simulate="intens",
                        radiusClust=NULL, clustering=5, verbose=TRUE) {

  #check if input arguments have correct values
  roundnumit <- round(numit)
  if ( roundnumit != numit ) {
    warning("numit must be a natural number. It is now rounded to: ", roundnumit)
    numit <- roundnumit
  }
  match.arg(simulate, choices=c("thinning", "intens", "clintens"))
  
  if (!is.null(hole)){
    winminus <- setminus.owin(A=ppdata$window, B=hole)
    ppdata <- ppp(x=ppdata$x, y=ppdata$y, window=winminus, marks=ppdata$marks)
  }

  if (!is.null(obsprobimage)){
    if (!is.null(hole)){
      holeim <- as.im(X=winminus, W=craterA$window, value=0)
      obsprobimage <- eval.im(obsprobimage + holeim)
    }
    sumim <- summary(obsprobimage)
    if(sumim$min <= 0){warning("obsprobimage contains value 0 or negative values")}
    if(simulate=="thinning"){warning("selected simulation scheme (thinning) does not take into account given obsprobimage")}
  }

  # weights which are used for intensity estimation later on (their reciprocal value)
  # if no obsprobimage is given, all marks (weights) equal 1
  if(!is.null(obsprobimage)){
    # Markierung mittels obsprobimage - hier noch nicht Kehrwert!
    obsprobs <- obsprobimage[ppdata]
  }else{obsprobs <- rep(1, times=ppdata$n)}
  invweight <- obsprobs

  if(min(invweight) <= 0){stop("obsprobs contains value 0 or negative values")}

  #here the intensity is being estimated
  if ( simulate == "intens" ) {

    origintens <- est_intens_weight(ppdata, covmatrix=covmatrix, weights=1/invweight)
    intensSim <- origintens$intensest
    intensSim$v <- (1/(1 - nxprob))*origintens$intensest$v

  } else if ( simulate == "clintens" ) {

    if( is.null(radiusClust) ) {
      radiusClust <- quantile(nndist(ppdata), p=0.7, type=8)
    }
    intensSim <- det_nsintens_restr(ppdata=ppdata, radius=radiusClust, weights=1/invweight)
  }

  result <- matrix(data=NA, nrow=0, ncol=6)

  numout <- 0
  i <- 1
  k <- 1
  alphastar <- cutoff


  while(i <= numit){

    if ( simulate == "thinning" ) {
      thinned <- thin(full=ppdata, nxprob=nxprob)
      if(!is.null(obsprobimage)){
        allobserved <- thinned$observed
        obsprobsobs <- obsprobimage[allobserved]
        indobs <- rbinom(n=allobserved$n, size=1, prob=obsprobsobs)
        observed <- allobserved[indobs==1]
      }else{observed <- thinned$observed}
      unobserved <- thinned$unobserved
    }
    if ( simulate == "intens" ) {
      thinned <- sim_intens(ppdata, intensSim, nxprob)
      if(!is.null(obsprobimage)){
        allobserved <- thinned$observed
        obsprobsobs <- obsprobimage[allobserved]
        indobs <- rbinom(n=allobserved$n, size=1, prob=obsprobsobs)
        observed <- allobserved[indobs==1]
      }else{observed <- thinned$observed}
      unobserved <- thinned$unobserved
    }
    if ( simulate == "clintens" ) {
      ppsim <- sim_nsprocess(ppdata=ppdata, intens=intensSim, radius=radiusClust,
                             clustering=clustering, thinning=nxprob)
      thinned <- thin(full=ppsim, nxprob=nxprob)
      if(!is.null(obsprobimage)){
        allobserved <- thinned$observed
        obsprobsobs <- obsprobimage[allobserved]
        indobs <- rbinom(n=allobserved$n, size=1, prob=obsprobsobs)
        observed <- allobserved[indobs==1]
      }else{observed <- thinned$observed}
      unobserved <- thinned$unobserved
    }


    if (is.null(intens)){
      if(!is.null(obsprobimage)){
        obsprobs <- obsprobimage[observed]
      }else{obsprobs <- rep(1, times=observed$n)}
      invweight <- obsprobs
      estim <- est_intens_weight(observed, covmatrix=covmatrix, weights=1/invweight)
      intens <- estim$intensest
      covmatrix <- estim$covmatrix
    }

    resultdetHRZ <- det_hrz_restr(ppdata=observed, type="intens", criterion="indirect",
                              cutoff=alphastar, hole=hole, integratehole=TRUE, obsprobimage=obsprobimage, intens=intens,
                              nxprob=nxprob, covmatrix=covmatrix, returnintens=FALSE)
    resultevalHRZ <- eval_hrz(hrz=resultdetHRZ$zone, unobspp=unobserved, obspp=observed)

    if(resultevalHRZ$numbermiss > 0){
      numout <- numout + 1
    }


    poutmin <- numout/numit
    poutmax <- (numout + numit - i)/numit

    resstep <- c(k, i, alphastar, numout, poutmin, poutmax)
    result <- rbind(result, resstep)

    if(poutmin > cutoff + tol){
      alphastar <- alphastar * i/(numit + 1)
      if(verbose){
        cat("Decrease alphastar to ", alphastar, " after ", i, " iterations with numout=", numout , "\n", sep="")
      }
      i <- 0
      numout <- 0
      k <- k + 1
    }

    if(poutmax < cutoff - tol){
      alphastar <- alphastar * (1 + (numit - i + 1)/numit)
      if(verbose){
        cat("Increase alphastar to ", alphastar, " after ", i, " iterations with numout=", numout , "\n", sep="")
      }
      i <- 0
      numout <- 0
      k <- k + 1
    }

    i <- i + 1

  }

  resultdf <- as.data.frame(result, row.names = as.character(1:dim(result)[1]))
  colnames(resultdf) <- c("k", "i", "alphastar", "numout", "poutmin", "poutmax")

  res <- list(alphastar=alphastar, course=resultdf)

  class(res) <- "bootcorr"
  return(res)

}
