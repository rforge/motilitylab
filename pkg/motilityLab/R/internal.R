.ulapply <- function( l, f, ... ) { 
  unlist(lapply(l,f, ...)) 
}


.minMaxOfTwoMatrices <- function(a, b, func) {
  #   print("a:")
  #   print(a)
  #   print("b:")
  #   print(b)
  res <- matrix(nrow=nrow(a), ncol=ncol(a))
  #   print(ncol(a))
  for (j in 1:ncol(a)) {
    #     print(j)
    res[1, j] <- min(a[1, j], b[1, j])
    res[2, j] <- max(a[2, j], b[2, j])
  }
  colnames(res) <- colnames(a)
  rownames(res) <- rownames(a)
  return(res)
}

.boundingBoxTrack <- function(track) {
  min.and.max.per.dim <- rbind(min = apply(track, 2, min), max = apply(track, 2, max))
  return(min.and.max.per.dim[, 2:ncol(min.and.max.per.dim)])
}

# TODO give min.segments a default value?
#' Compute a Measure in a Segmentwise Fashion
#' 
#' For each number of segments the mean is computed over the measure applied to
#' all the track's subtracks consisting of as many single segments.  
#' 
#' @param track the (super-) track whose subtracks shall be regarded, given in 
#' a data frame with columns \eqn{t}, \eqn{x} and possibly \eqn{y}, and \eqn{z}.
#' @param measure the measure that is to be computed for all subtracks
#' @param min.segments the minimum number of segments of which the
#' regarded subtracks shall consist. For example when computing a measure that 
#' compares segments, e.g. the \code{\link{overallAngle}}, this value should 
#' be set to (at least) 2, otherwise it should be (at least) 1. 
#' @param only.for.i either a numeric \eqn{i} indicating that only subtracks of
#' \eqn{i} segments are to be regarded or \code{NULL} to express that all 
#' subtracks that constist of at least \code{min.segments} and at most \eqn{n-1}
#' segments, where \eqn{n} is the number of points in the \code{track}, are to 
#' be considered.
#' 
#' @return With \eqn{n} being the number of points in the input track, a vector
#' of length \deqn{n - min.segments} ist returned, containing in 
#' the \eqn{i}th position the mean value of the measure for subtracks 
#' consisting of \eqn{i} segments.
#' 
#' @details For each number of segments at least as large as a lower bound given in the 
#' third parameter, this function computes the mean of the measure applied on
#' every subtrack of at least so many single segments.
.computeSegmentwiseMeans <- function(track, measure, min.segments=1, only.for.i=NULL, ...) { 
  if (nrow(track) <= min.segments) {
    return(NA)
  }
  means <- c()
  debug.nr.tracks <- 0
  if (is.null(only.for.i)) {
    for (i in (min.segments):(nrow(track) - 1)) {
      subtracks <- computeSubtracksOfISegments(track, i, ...)
      val <- sapply(subtracks, measure)
      means[i - min.segments + 1] <- mean(val)#, na.omit=TRUE)    # na.rm produces NaN if all values are NA
    }
  } else {
    if ((only.for.i >= min.segments) && (only.for.i <= (nrow(track) - 1))) {
      subtracks <- computeSubtracksOfISegments(track, only.for.i)
      val <- sapply(subtracks, measure)
      means <- mean(val)#, na.omit=TRUE)
    } else {
      return(NA)
    }
  }
  return(means)
}


# TODO wenn man die fkt turning angle bekommt automatisch mean.seg auf 2 setzen
# TODO segmentwise doesn't really fit, maybe rather "segementLengthWise"?
#' Segmentwise Version of a Measure
#' 
#' Returns a function which computes a the input measure segmentwise 
#' (see \code{\link{.computeSegmentwiseMeans}}) for a  given track.
.segmentwiseMeans <- function(measure, ...) {
  return(function(track) {
    .computeSegmentwiseMeans(track, measure, ...) #TODO DOk: hier fehlt noch die min.length.of.segments!!
  })
}



