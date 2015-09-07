#' Track Length
#' 
#' Estimates the length of the track as the sum of the lengths of its
#' linear segments. (Note that this is an underestimation).
#'
#' @param track The track whose length is to be computed.
#' @return The length.
#' @details The track length is computed by summing the Euclidean distances 
#' between consecutive positions.
trackLength <- function(track) {
	if (nrow(track) > 2) { 
		dif <- apply(track[,-1], 2, diff)
		return(sum(sqrt(apply(dif^2, 1, sum))))
	} else if (nrow(track) == 2) {
		# this case is necessary because if dimension dropping by 'apply'
		return(sqrt(sum((track[2,-1] - track[1,-1])^2)))
	} else if (nrow(track) == 1) {
		return(0)
	} else {
		return(NA)
	}
}


#' Mean Track Speed
#' 
#' Estimate the mean speed of a track by dividing its 
#' length by its duration.
#' 
#' @param track a track.
#' @return an estimate of the track's average speed.
#' @details The track's mean speed is computed by dividing its 
#' \code{\link{trackLength}} by its \code{\link{duration}}.
#' @seealso \code{\link{trackLength}}, \code{\link{duration}}
speed <- function(track) {
  trackLength(track) / duration(track)
}

#' Track Duration
#' 
#' Returns the duration of a track.
#' 
#' @param track a track.
#' 
#' @details Computes the track's duration by subtracting the 
#' time of the its first point from the time of its last.
#' 
#' @return Returns a tracks duration.
duration <- function(track) {
  dur <- track[nrow(track), 1] - track[1,1]
  dur <- unname(dur)
  return(dur)
}

#' Track Displacement
#' 
#' Computes the distance between the track's first and final positions.
#'
#' @param track the track whose displacement is to be computed.
#' @param limits Vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.

#' @return The track's displacement will be returned.
#' @details Computes the Euclidean distance between the track's start and end points.
#'
displacement <- function(track, limits=c(1,nrow(track))) {
  sqrt(sum((track[limits[2], -1] - track[limits[1], -1])^2))
}

#' Vector Between Track Endpoints
#' 
#' Computes the coordinates of the vector between that tracks' start and end points.
#' 
#' @param track the track whose displacement vector is to be computed.
#' 
#' @details The coordinates of the tracks's starting point are subtracted 
#' from its end point's coordinates.
#' 
#' @return The difference between the coordinates of the track's end point and 
#' those of its starting point is returned.
displacementVector <- function(track) {
  ret <- track[nrow(track), 2:ncol(track)] - track[1, 2:ncol(track)]
  rownames(ret) <- NULL  
  return(as.vector(ret))
}

#' Maximal Displacement
#' 
#' Computes the maximal distance of any position on the track from the 
#' starting position.
#' 
#' @param track the track whose maximal displacement is to be computed.
#' @param limits Vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' @return The maximal displacement.
#' @details Computes the maximum over all pairwise Euclidean distances between 
#' the track's starting point and another point on the track.
maxDisplacement <- function(track, limits=c(1,nrow(track))) {
	sqrt(max(rowSums(sweep(track[seq(limits[1],limits[2]),-1],2,track[limits[1],-1])^2)))
}

#' Square Displacement
#' 
#' Computes the squared distance between the track's first and final positions.
#'
#' @param track the track whose displacement is to be computed.
#' @param limits Vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' 
#' @return The track's square displacement will be returned.
#' 
#' @details Computes the squared Euclidean distance between the track's start 
#' and end point, by not extrakting the square root of the squared distances on 
#' the single dimensions.
#'
squareDisplacement <- function(track, limits=c(1,nrow(track))) {
  sum((track[limits[2], -1] - track[limits[1], -1])^2)
}

#' Track Aspericity
#' 
#' Computes the asphericity of the track's points.
#' 
#' @param track the track whose asphericity is to be computed. Note that the
#' track must be in a two- or threedimensional space, thus the data frames 
#' containing the track's points must contain an \eqn{x} and \eqn{y} value. 
#' @return The track's asphericity, a vaule between \eqn{0} and \eqn{1}, will 
#' be returned for two- or threedimensional tracks, otherwise \code{NA}.
#' @details Computes the asphericity of a track, i.e. the scatter plot of 
#' its points, by considering the length of its principal components, 
#' see References. The value will be from the interval \eqn{[0,1]}, 
#' with high values indicating track straightness. The asphericity is a 
#' measure for the track's straightness. Note that the asphericity is only 
#' defined for tracks in two or three dimensions.
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
asphericity <- function(track) {
	dim <- ncol(track) - 1
  if (dim == 1) {
    return(NA)
  }
  if (nrow(track)==1) {
    return(1)
  }
	eigen.values <- eigen(cov(track[,2:ncol(track)]))$values
	rav <- mean(eigen.values)
	res <- sum((eigen.values - rav )^2 / dim / (dim - 1) / rav^2)
#   if (is.nan(res)) print(track)
  return(res)
}


#' Track Straightness
#' 
#' Computes a track's straightness, i.e., the track's displacement divided by 
#' its length.
#'
#' @param track the track whose straightness is to be computed.
#' @return The track's straightness, a value between \eqn{0} and \eqn{1}, will 
#' be returned. If the track has length \eqn{0}, then NaN is returned.
#' @details Computes the straightness of a track, i.e. the ratio of the 
#' track's displacement and length. The value will be from the interval 
#' \eqn{[0,1]}, or 1 if the track has length \eqn{0}.
#' @seealso \code{\link{displacement}}, \code{\link{trackLength}}
straightness <- function(track) {
  n <- trackLength(track)
  if (n > 0) {
    return(displacement(track) / n)
  } else {
    return(1)
  }
}

#' Track Displacement Ratio
#' 
#' Computes a track's displacement ratio, i.d. the track's (current) 
#' displacement divided by its maximal displacement.
#'
#' @param track the track whose displacement ratio is to be computed.
#' @return The track's displacement ratio, a vaule between \eqn{0} and \eqn{1}, 
#' will be returned. If the track has a maximal displacement of \eqn{0}, then 
#' NaN is returned.
#' @details Computes the displacement ratio of a track, i.e. the ratio of the 
#' track's (current) displacement and its maximum displacement. The value will
#' be from the interval \eqn{[0,1]}, or 1 if the track has a maximal 
#' displacement of 0.
#' @seealso \code{\link{displacement}}, \code{\link{maxDisplacement}}
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
displacementRatio <- function(track) {
  dmax <- maxDisplacement(track)
  if (dmax > 0) {
    return(displacement(track) / dmax)
  } else {
    return(1)
  }
}

#' Track Outreach Ratio
#' 
#' Computes a track's maximal displacement divided by its length.
#'
#' @param track the track whose outreach ratio is to be computed.
#' @return The track's outreach ratio, a value between \eqn{0} and \eqn{1}, 
#' will be returned. If the track has length \eqn{0}, then NaN is returned.
#' @details Computes the outreach ratio of a track, i.e. the ratio of the 
#' track's maximal displacement and length. The value will be from the interval 
#' \eqn{[0,1]}, or 1 if the track has length \eqn{0}.
#' @seealso \code{\link{maxDisplacement}}, \code{\link{trackLength}}
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
outreachRatio <- function(track) {
  trackLength <- trackLength(track) 
  if (trackLength > 0) {
    return(maxDisplacement(track) / trackLength)
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
    return( acos(sum(a * b) / (sqrt(sum(a^2) * sum(b^2)))) )
  }
}

#' Dot Product Between Steps at Track Endpoints
#'
#' Computes the dot product between the first and last steps of the given track.
#' @param x the input track.
#' @param limits vector giving the first and last row of the track. Can be used to avoid
#' extracting subtracks, which is exploited e.g. by 
#' \code{\link{aggregate.tracks}}.
#' @return the dot product. 
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
#' @param track the track whose mean turning angle is to be computed.
#' @details Computes the angle between each two subsequent segments in the 
#' track and returns their mean. Angles are metered symmetrically, thus 
#' yielding (degree) values between \eqn{0} and \eqn{180}. (Both a \eqn{90} 
#' degrees left and right turn yield the value 90.)
#' @return The tracks mean turning angle will be returned in degrees. 
meanTurningAngle <- function(track) {
	mean(sapply(subtracks(track, 2), overallAngle), na.rm=TRUE)
}


#' Track Hurst Exponents
#'
#' Computes the Hurst exponent of each of the track's dimensions.
#' 
#' @param track the track whose Hurst exponent is to be computed.
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
hurstExponent <- function(track) {
	if( !requireNamespace("pracma",quietly=TRUE) ){
		stop("This function requires the 'pracma' package.")
	}
  if (nrow(track) < 6) {
    return(NA)
  }
  if ((nrow(track)%%2)==1) {
    track <- track[1:(nrow(track)-1), ]
  }
  he <- function(x) {
    pracma::hurstexp(x, display=FALSE)$Hs
  }
  ret <- c()
  ret <- apply(as.matrix(track[,2:ncol(track)]), 2, he)
  return(ret)
}


#' Track Fractal Dimension
#'
#' Computes the fractal dimension of a track using all track dimensions by the 
#' box-count method.
#'
#' @param track the track for which the fractal dimension is to be calculated.
#' @details The fractal dimension is estimated using the function from
#' \code{\link[fractaldim]{fd.estim.boxcount}} from the 
#' \code{fractaldim} package. For Brownian motion, fractal dimension and Hurst exponent
#' (see \link{hurstExponent}) 
#' are related by the formula H=2-D. For non-Brownian motion, however, this relationship
#' need not hold. While the hurst exponent takes a global approach to the 
#' track's properties, fractal dimension is a local approach to the track's properties
#' (Gneiting and Schlather, 2004).
#'
#' @return the number indicating the track's fractal dimension.
#'
#' @seealso \code{\link{hurstExponent}}
#'
#' @references 
#' Tillmann Gneiting and Martin Schlather (2004), Stochastic Models That Separate Fractal
#' Dimension and the Hurst Effect. \emph{SIAM Review} \bold{46}(2), 269--282. 
#' doi:10.1137/S0036144501394387
#'
fractalDimension <- function(track){
  if( !requireNamespace("fractaldim",quietly=TRUE) ){
    stop("This function requires the 'fractaldim' package.")
  }
  return(fractaldim::fd.estim.boxcount(track[,-1])$fd)
}
