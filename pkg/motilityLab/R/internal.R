.ulapply <- function( l, f, ... ) { 
  unlist(lapply(l,f, ...), use.names=FALSE) 
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

.subtrackIndices <- function( track, i, overlap ){
  if (nrow(track) <= i) {
    return(NA)
  } 
  if (overlap > i-1) {
    warning("Overlap exceeds segment length")
    overlap <- i-1
  }
  x <- seq(1, (nrow(track)-i),i-overlap)
  cbind( first=x, last=x+i )
}


.computeSegmentwiseMeans <- function(track, measure, min.segments=1, only.for.i=NULL, ...) { 
  if (nrow(track) <= min.segments) {
    return(NA)
  }
  means <- c()
  debug.nr.tracks <- 0
  if (is.null(only.for.i)) {
    for (i in (min.segments):(nrow(track) - 1)) {
      subtracks <- .computeSubtracksOfISegments(track, i, ...)
      val <- sapply(subtracks, measure)
      means[i - min.segments + 1] <- mean(val)#, na.omit=TRUE)    # na.rm produces NaN if all values are NA
    }
  } else {
    if ((only.for.i >= min.segments) && (only.for.i <= (nrow(track) - 1))) {
      subtracks <- .computeSubtracksOfISegments(track, only.for.i)
      val <- sapply(subtracks, measure)
      means <- mean(val)#, na.omit=TRUE)
    } else {
      return(NA)
    }
  }
  return(means)
}

.computeSubtracksOfISegments <- function(track, i, overlap=i-1) {
  if (nrow(track) <= i) {
    return( NULL )
  } 
  if (overlap > i-1) {
    warning("Overlap exceeds segment length")
    overlap <- i-1
  }
  l <- list()
  for (j in seq(1, (nrow(track)-i),i-overlap)) {
    l[[as.character(j)]] <- track[j:(j+i), ,drop=FALSE]
  }    
  structure(l, class="tracks")
}


.segmentwiseMeans <- function(measure, ...) {
  return(function(track) {
    .computeSegmentwiseMeans(track, measure, ...) #TODO DOk: hier fehlt noch die min.length.of.segments!!
  })
}



