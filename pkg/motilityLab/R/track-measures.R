#' Track Length
#' 
#' Estimates the length of the track as the sum of the lengths of its
#' linear segments. 
#'
#' @param x The track whose length is to be computed.
#' @return The length.
#' @details This amounts to linar interpoplation, which really gives a lower
#' bound on the track length and is therefore usually an underestimation.
trackLength <- function(x) {
	if (nrow(x) > 2) { 
		dif <- apply(x[,-1], 2, diff)
		return(sum(sqrt(apply(dif^2, 1, sum))))
	} else if (nrow(x) == 2) {
		# this case is necessary because if dimension dropping by 'apply'
		return(sqrt(sum((x[2,-1] - x[1,-1])^2)))
	} else if (nrow(x) == 1) {
		return(0)
	} else {
		return(NA)
	}
}


#' Mean Track Speed
#' 
#' Estimate the mean instantaneous speed of a track by dividing its 
#' \code{\link{trackLength}} by its \code{\link{duration}}.
#' 
#' @param x a track.
#' @return An estimate of the track's average speed, a nonnegative number.
speed <- function(x) {
  trackLength(x) / duration(x)
}

#' Track Duration
#' 
#' Returns the duration of a track, i.e., the time elapsed between its first and last 
#' position.
#' 
#' @param x a track.
#' @return The track's duration, a nonnegative number.
duration <- function(x) {
  dur <- x[nrow(x), 1] - x[1,1]
  dur <- unname(dur)
  return(dur)
}

#' Track Displacement
#' 
#' Computes the Euclidean distance between the track's start and end points.
#'
#' @param x the track whose displacement is to be computed.
#' @param limits Vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.

#' @return The track's displacement, a nonnegative number.
#'
displacement <- function(x, limits=c(1,nrow(x))) {
  sqrt(sum((x[limits[2], -1] - x[limits[1], -1])^2))
}

#' Vector Between Track Endpoints
#' 
#' Computes the vector between the track endpoints.
#' 
#' @param x the track whose displacement vector is to be computed.
#' @return The displacement vector, which has as many dimensions as the input tracks.
displacementVector <- function(x) {
  ret <- x[nrow(x),-1] - x[1,-1]
  rownames(ret) <- NULL  
  return(as.vector(ret))
}

#' Maximal Displacement
#' 
#' Computes the maximal distance of any position on the track from the 
#' starting position.
#' 
#' @param x the track whose maximal displacement is to be computed.
#' @param limits vector giving the first and last row of the track. Can be used to avoid
#' computing subtracks explicitly, which is exploited by 
#' \code{\link{aggregate.tracks}}.
#' @return The maximal displacement, a nonnegative number.
maxDisplacement <- function(x, limits=c(1,nrow(x))) {
	sqrt(max(rowSums(sweep(x[seq(limits[1],limits[2]),-1],2,x[limits[1],-1])^2)))
}

#' Square Displacement
#' 
#' Computes the squared distance between the track's first and final positions.
#'
#' @param x the track whose displacement is to be computed.
#' @param limits Vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' 
#' @return The track's square displacement will be returned.
#' 
squareDisplacement <- function(x, limits=c(1,nrow(x))) {
  sum((x[limits[2], -1] - x[limits[1], -1])^2)
}

#' Track Aspericity
#' 
#' Computes the asphericity of the set of positions on the track
#' via the length of its principal components
#' (Mokhtari et al, 2013).  Asphericity is a 
#' measure for the track's straightness: it is a number between 0 and 1, with 
#' higher values indicating straighter tracks. 
#' In that sense it is similar to \code{\link{straightness}}, however, the asphericity is 
#' not affected by back-and-forth motion of the object. 
#' 
#' @param x input track, which must have two or more spatial dimensions.
#' @param limits vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' @return The track's asphericity, a value between \eqn{0} and \eqn{1}, is 
#' returned for two- or higher dimensional tracks. Otherwise the return value is \code{NA}.
#' @details We define the asphericity of every track with two or fewer positions to be 1. 
#' This is true even in the case where two positions are identical, which is motivated 
#' by considering
#' that situation as the limit of a process in which two points move towards each other.
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
#'
#' @seealso \code{\link{straightness}}
#'
asphericity <- function(x,limits=c(1,nrow(x))) {
  dim <- ncol(x) - 1
  if (dim == 1) {
    return(NA)
  }
  if (limits[2]-limits[1]<3) {
    return(1)
  }
  eigen.values <- eigen(cov(x[limits[1]:limits[2],-1]))$values
  rav <- mean(eigen.values)
  res <- sum((eigen.values - rav )^2 / dim / (dim - 1) / rav^2)
  return(res)
}


#' Track Straightness
#' 
#' Computes a track's straightness, i.e., the track's \code{\link{displacement}} 
#' divided by its \code{\link{trackLength}}.
#'
#' @param x the track whose straightness is to be computed.
#' @return The track's straightness, a value between \eqn{0} and \eqn{1}. 
#' If the track has \eqn{0}, then 1 is returned.
#' @seealso \code{\link{asphericity}} for a different measure of straightness.
straightness <- function(x) {
  l <- trackLength(x)
  if (l > 0) {
    return(displacement(x) / l)
  } else {
    return(1)
  }
}

#' Track Displacement Ratio
#' 
#' Computes a track's displacement ratio, i.d. the track's (current) 
#' displacement divided by its maximal displacement.
#'
#' @param x the track whose displacement ratio is to be computed.
#' @return The track's displacement ratio, a value between \eqn{0} and \eqn{1}, 
#' will be returned. If the track has a maximal displacement of \eqn{0}, then 
#' 1 is returned.
#' @seealso \code{\link{displacement}}, \code{\link{maxDisplacement}}
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
displacementRatio <- function(x) {
  dmax <- maxDisplacement(x)
  if (dmax > 0) {
    return(displacement(x) / dmax)
  } else {
    return(1)
  }
}

#' Track Outreach Ratio
#' 
#' Computes a track's maximal displacement divided by its length.
#' The value will be from the interval \eqn{[0,1]}, or 1 if the track has length 0.
#'
#' @param x the track whose outreach ratio is to be computed.
#' @return The track's outreach ratio, a value between \eqn{0} and \eqn{1}, 
#' will be returned. If the track has length \eqn{0}, then NaN is returned.
#' @details Computes the outreach ratio of a track, i.e. the ratio of the 
#' track's maximal displacement and length. 
#' @seealso \code{\link{maxDisplacement}}, \code{\link{trackLength}}
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
outreachRatio <- function(x) {
  l <- trackLength(x) 
  if (l > 0) {
    return(maxDisplacement(x) / l)
  } else {
    return(1)
  }
}

#' Angle Between Steps at Track Endpoints
#' 
#' Computes the angle between the first and the last segment of the given track.
#' @param x the track whose overall turning angle is to be computed.
#' @param limits vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' @details Computes the angle between the vectors representing the track's 
#' first and last segment, respectively, i.e. the overall turning angle.
#' Angles are measured symmetrically, thus the return values range from 0 to 180 degrees.
#' For instance, both a 90 degrees left and right turn yield the value 90.
#'
#' @return The track's overall turning angle in radians.
#' @examples
#' ## show a turning angle plot with error bars for the T cell data.
#' with( (aggregate(BCells,overallDot,FUN="mean.se",na.rm=TRUE)),{
#'   plot( mean ~ i, xlab="time step", 
#'   	ylab="turning angle (rad)", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
overallAngle <- function(x, limits=c(1,nrow(x))) {
  if (limits[2]-limits[1] < 2) {
    return(0)
  } else {
    a <- diff(x[limits[1]:(limits[1]+1),-1])
    b <- diff(x[(limits[2]-1):limits[2],-1])
    a <- a/sqrt(sum(a^2))
    b <- b/sqrt(sum(b^2))
    return( acos(sum(a * b)) )
  }
}

#' Dot Product Between Steps at Track Endpoints
#'
#' Computes the dot product between the first and last steps of the given track.
#' @param x the input track.
#' @param limits vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' @return The dot product. 
#' @examples
#' ## compute and plot the autocovariance function for the T cell data (assuming isotropy)
#' with( (aggregate(BCells,overallDot,FUN="mean.se",na.rm=TRUE)),{
#'   plot( cos(mean) ~ i, xlab="time", 
#'   ylab="correlation of orientation", type="l" )
#'   segments( i, cos(lower), y1=cos(upper) )
#' } )
overallDot <- function(x, limits=c(1,nrow(x))) {
  a <- diff(x[limits[1]:(limits[1]+1),-1])
  b <- diff(x[(limits[2]-1):limits[2],-1])
  return( sum(a * b) )
}


#' Track Mean Turning Angle
#' 
#' Computes the mean of all angles between two subsequent segments of the given 
#' track.
#' @param x the track whose mean turning angle is to be computed.
#' @details Computes the angle between each two subsequent segments in the 
#' track and returns their mean. Angles are metered symmetrically, thus 
#' yielding (degree) values between \eqn{0} and \eqn{180}. (Both a \eqn{90} 
#' degrees left and right turn yield the value 90.)
#' @return The tracks mean turning angle will be returned in degrees. 
meanTurningAngle <- function(x) {
	mean(sapply(subtracks(x, 2), overallAngle), na.rm=TRUE)
}


#' Track Hurst Exponents
#'
#' Computes the Hurst exponent of each of the track's dimensions.
#' 
#' @param x the track whose Hurst exponent is to be computed.
#' @details If the track includes at least six and an even number of points,
#' for each of the track's dimensions, the empirical Hurst exponent is 
#' determined using the function \code{\link[pracma]{hurstexp}} from the 
#' \code{pracma} package and returned as part of a vector.
#' If the track consists of less than six points, NA is returned, if it 
#' comprises an uneven number of points, the last one is discarded.
#' @return A vector with the Hurst exponent for each of the tracks dimension 
#' is returned. If the track has an uneven number of points, the last one is 
#' discarded, if it has less than six points, NA will be returned.
#'
#' @seealso \code{\link{fractalDimension}}
#'
hurstExponent <- function(x) {
	if( !requireNamespace("pracma",quietly=TRUE) ){
		stop("This function requires the 'pracma' package.")
	}
  if (nrow(x) < 6) {
    return(NA)
  }
  if ((nrow(x)%%2)==1) {
    x <- x[1:(nrow(x)-1), ]
  }
  he <- function(x) {
    pracma::hurstexp(x, display=FALSE)$Hs
  }
  ret <- c()
  ret <- apply(as.matrix(x[,-1]), 2, he)
  return(ret)
}


#' Track Fractal Dimension
#'
#' Computes the fractal dimension of a track using all track dimensions by the 
#' box-count method.
#'
#' @param x the track for which the fractal dimension is to be estimated.
#' @details The fractal dimension is estimated using the function from
#' \code{\link[fractaldim]{fd.estim.boxcount}} from the 
#' \code{fractaldim} package. For Brownian motion, fractal dimension and Hurst exponent
#' (see \link{hurstExponent}) 
#' are related by the formula H=2-D. For non-Brownian motion, however, this relationship
#' need not hold. While the hurst exponent takes a global approach to the 
#' track's properties, fractal dimension is a local approach to the track's properties
#' (Gneiting and Schlather, 2004).
#'
#' @return The estimate of the track's fractal dimension.
#'
#' @seealso \code{\link{hurstExponent}}
#'
#' @references 
#' Tillmann Gneiting and Martin Schlather (2004), Stochastic Models That Separate Fractal
#' Dimension and the Hurst Effect. \emph{SIAM Review} \bold{46}(2), 269--282. 
#' doi:10.1137/S0036144501394387
#'
fractalDimension <- function(x){
  if( !requireNamespace("fractaldim",quietly=TRUE) ){
    stop("This function requires the 'fractaldim' package.")
  }
  return(fractaldim::fd.estim.boxcount(x[,-1])$fd)
}
