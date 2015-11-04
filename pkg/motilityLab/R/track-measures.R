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
#' @param from index of first row of the track. 
#' @param to index of last row of the track. 
#' Both `from` and `to` can be vectors and are useful to avoid extracting subtracks, 
#' which is exploited by \code{\link{aggregate.tracks}}.
#' 
#' @return The track's displacement, a nonnegative number.
#'
displacement <- function( x, from=1, to=nrow(x) ) {
  sqrt(.rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
  	length(from),ncol(x)-1))
}

#' Vector Between Track Endpoints
#' 
#' Computes the vector between the track endpoints.
#' 
#' @param x the input track.
#'
displacementVector <- function(x) {
  ret <- x[nrow(x),-1] - x[1,-1]
  rownames(ret) <- NULL  
  return(as.vector(ret))
}

#' Maximal Displacement
#' 
#' Computes the maximal Euclidean distance of any position on the track from the 
#' first position.
#' 
#' @param x the input track.
#' @return The maximal displacement, a nonnegative number.
maxDisplacement <- function(x) {
	limits <- c(1,nrow(x))
	sqrt(max(rowSums(sweep(x[seq(limits[1],limits[2]),-1],2,x[limits[1],-1])^2)))
}

#' Squared Displacement
#' 
#' Computes the squared distance between the track's first and final positions.
#'
#' @param x the input track.
#' @param from index of first row of the track. 
#' @param to index of last row of the track. 
#' Both `from` and `to` can be vectors and are useful to avoid extracting subtracks, 
#' which is exploited by \code{\link{aggregate.tracks}}.
#' 
#' @return The track's square displacement will be returned.
#' 
squareDisplacement <- function(x, from=1, to=nrow(x)) {
  .rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
  	length(from),ncol(x)-1)
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
asphericity <- function(x) {
  dim <- ncol(x) - 1
  if (dim == 1) {
    return(NA)
  }
  limits <- c(1,nrow(x))
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
#' @param from index of first row of the track. 
#' @param to index of last row of the track. 
#' Both `from` and `to` can be vectors and are useful to avoid extracting subtracks, 
#' which is exploited by \code{\link{aggregate.tracks}}.
#' @param xdiff row differences of x, should be pre-computed if this function is calleed
#' repeatedly for subtracks of a long track.
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
overallAngle <- function(x, from=1, to=nrow(x), xdiff=diff(x)) {
	r <- rep(0, length(from))
	ft <- from<(to-1)
	a <- xdiff[from[ft],-1,drop=FALSE]
	b <- xdiff[to[ft]-1,-1,drop=FALSE]
	a <- a/sqrt(.rowSums(a^2, nrow(a), ncol(a)))
	b <- b/sqrt(.rowSums(b^2, nrow(b), ncol(b)))
	r[ft] <- acos(.rowSums(a * b, nrow(a), ncol(a)))
	r
}

#' Dot Product Between Steps at Track Endpoints
#'
#' Computes the dot product between the first and last steps of the given track.
#' @param x the input track. Must have at least 2 rows or an error is thrown.
#' @param from index of first row of the track. 
#' @param to index of last row of the track. 
#' Both `from` and `to` can be vectors and are useful to avoid extracting subtracks, 
#' which is exploited by \code{\link{aggregate.tracks}}.
#' @param xdiff row differences of x, should be pre-computed if this function is calleed
#' repeatedly for subtracks of a long track.
#' @return The dot product. 
#' @examples
#' ## compute and plot the autocovariance function for the T cell data 
#' ## (assuming isotropy)
#' with( (aggregate(BCells,overallDot,FUN="mean.se")),{
#'   plot( mean ~ i, xlab="time", 
#'   ylab="autocovariance", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
overallDot <- function(x, from=1, to=nrow(x), xdiff=diff(x)) {
	r <- rep(NaN, length(from))
	ft <- from<to
	a <- xdiff[from[ft],-1,drop=FALSE]
	b <- xdiff[to[ft]-1,-1,drop=FALSE]
	r[ft] <- .rowSums(a * b, nrow(a), ncol(a))
	r
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


#' Track Hurst Exponent
#'
#' Computes the empirical Hurst exponent of the track.
#' 
#' @param x the input track.
#' @return The corrected empircal Hurst exponent as estimated by 
#' \code{\link[pracma]{hurstexp}}
#' is returned. If the track has less than two 
#' positions, NA will be returned.
#'
#' @seealso A related measure is the \code{\link{fractalDimension}}.
#'
hurstExponent <- function(x) {
	if( !requireNamespace("pracma",quietly=TRUE) ){
		stop("This function requires the 'pracma' package.")
	}
  if (nrow(x) < 2) {
    return(NA)
  }
  return(pracma::hurstexp(x[,-1], display=FALSE)$Hal)
}


#' Track Fractal Dimension
#'
#' Computes the fractal dimension of a track using all track dimensions by the 
#' box-count method.
#'
#' @param x the input track.
#' @details The fractal dimension is estimated using the function
#' \code{\link[fractaldim]{fd.estim.boxcount}} from the 
#' \code{fractaldim} package. For self-affine processes in n dimensions, 
#' fractal dimension and Hurst exponent
#' (see \code{\link{hurstExponent}}) 
#' are related by the formula H=n+1-D. 
#' For non-Brownian motion, however, this relationship
#' need not hold. While the Hurst exponent takes a global approach to the 
#' track's properties, fractal dimension is a local approach to the track's properties
#' (Gneiting and Schlather, 2004).
#'
#' @return The estimate of the track's fractal dimension.
#'
#' @seealso A related measure is the \code{\link{hurstExponent}}.
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
