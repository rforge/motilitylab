# TODO Johannes: Stimmt die Doku?
#' Convert from Tracks to Data Frame
#' 
#' Converts a \code{tracks} object into a data frame.
#' 
#' @param x the  \code{tracks} object to be coerced to a data frame.
#' 
#' @param row.names ‘NULL’ or a character vector giving the row names for the
#'  data frame.  Missing values are not allowed.
#'
#' @param optional logical. Column names are always assigned to the resulting
#'  data frame regardless of the setting of this option.
#'
#' @details Returns one data frame containing all indiviual 
#' tracks from the input with a prepended column "id" containing
#' each track's name in the list of the \emph{tracks} object.
#' 
#' @return Returns one data frame containing all data of the indiviual 
#' tracks from the input with a prepended column named "id" which contains 
#' each track's name in the list of the \emph{tracks} object.
as.data.frame.tracks <- function(x, row.names = NULL, optional = FALSE, ...) {
	ids <- rep(names(x), lapply(x,nrow))
	r <- cbind(id=ids, Reduce(rbind, x))
	if( !is.null(row.names) && length(row.names)==nrow(r) ){
		rownames(r) <- row.names
	}
	return(r)
}

#' Convert to tracks
#' 
#' Coerces \code{x} to a \code{tracks} object.
#' 
#' @param x the object to be coerced.
#'
#' @details the S3 Method corresponding to \code{x}'s type (usually a list) is called.
#'  
#' @return the coerced object of S3 class \code{tracks}.
as.tracks <- function(x,...) 
  UseMethod("as.tracks")


#' Convert from Tracks to List
#' 
#' Coerces a \code{tracks} object to a list
#' 
#' @param x the \code{tracks} object to be coerced to a list.
#' 
#' @details the input \code{tracks} object is coerced to a list.
#' 
#' @return \code{x} is returned as a list (containing the
#' data frames which represent the tracks)
as.list.tracks <- function(x, ...) 
  structure( x, class="list" )


#' Convert from List to Tracks
#' 
#' Converts a List into a \emph{Tracks} Object
#' 
#' @param x the list that is to be converted into a \emph{tracks} object.
#' 
#' @details The input list is assigned the S3 class \code{tracks}.
#' 
#' @return the coerced object of S3 class \code{tracks} (i.e., a list 
#'  containing a data frame for each track)
as.tracks.list <- function(x,...) 
  structure(x, class="tracks")


## TODO remove for loop
#' Get a Subset of Tracks
#' 
#' Returns a \emph{tracks} object containing only the specified elements of a given
#' \emph{tracks} object.
#' 
#' @param tracks the tracks object from which certain tracks are to be chosen.
#' @param ids a vector containing the ids of the tracks that shall be chosen.
#' 
#' @details The elements given in \code{ids} are chosen from the list that represents
#' the \emph{tracks} object and returned, likewise as a \emph{tracks} object.
#' 
#' @return A \emph{tracks} object containing just the elements from the given 
#' \emph{tracks} object which are specified in \code{ids}.
getTracks <- function(tracks, ids) {
	r <- list()
	for (i in ids) {
		r[[i]] <- tracks[[i]]
	}
	structure(r, class="tracks")
}


#' Sort each Track's Positions by Time
#' 
#' Sorts the positions in each track in a \emph{tracks} object by time.
#' 
#' @param x the \emph{tracks} object whose tracks are to be sorted by time.
#' 
#' @param decreasing logical.  Should the sort be increasing or decreasing? 
#'  Provided only for consistency with the generic sort method. The positions in
#'  each track should be sorted in increasing time order.
#'
#' @details Sorts the positions of each track (represented as a data frame) in the 
#' \emph{tracks} object by time (given in the column t). 
#' 
#' @return A \emph{tracks object} that contains the tracks from the input object 
#' sorted by time is returned.
sort.tracks <- function(x, decreasing=FALSE, ...) {
	as.tracks(lapply(x, function(d) d[order(d$t),decreasing=decreasing] ))
}


#' Concatenate \emph{Tracks} Objects
#' 
#' Concatenates several \emph{tracks} objects.
#' 
#' @param ... the \emph{tracks} objects that are to be concatenated.

#' @details Joins the tracks (given as data frames) in the given \emph{tracks} 
#' objects into one tracks object, preserving the names of list elements and 
#' the ordering of the input objects as well as of the individual tracks in 
#' them.
#' 
#' @return A \emph{tracks} object containing all the tracks given in any of 
#' the input \emph{tracks} objects.
c.tracks <- function(...) {
	args <- lapply(list(...), as.list)
	as.tracks(do.call(c, args))
}


# TODO: What is pos?
#' Data Input form a CSV File
#' 
#' Reads cell tracks from a CSV file which contains a number of tracks, each 
#' with a distinct value in the first field 'id', and with at least 3 points 
#' per track, represented by their pos (?), time stamp \eqn{t} and 
#' coordinate(s), \eqn{x}, \eqn{x} and \eqn{y} or \eqn{x}, \eqn{y} and \eqn{z}.
#' 
#' @param file the name of the file which the data are to be read from, a 
#' readable text-mode connection or a complete URL
#' (see \code{\link[utils]{read.table}}).
#' 
#' @param ... Further arguments to be passed to \code{read.csv}.
#' 
#' @return An object of class \emph{tracks} is returned, which is a list of 
#' data frames, each containing the information on one track. The data frames 
#' have a column \eqn{t}, followed by one column for each of the input track's 
#' coordinates.
#' 
#' @details The input file's first four fields are interpreted as \eqn{id}, 
#' \eqn{pos}, \eqn{t} and \eqn{x}, respectively, and, if available, the fifth
#' as \eqn{y} and the sixth as \eqn{z}. The returned object has the class 
#' \emph{tracks}, which is a list of data frames representing the single 
#' tracks and having columns \eqn{t} and \eqn{x}, plus \eqn{y} and \eqn{z}, if
#' necessary. The tracks' ids are retained in their position in the list, while
#' the field \eqn{pos} will be unmaintained.
read.tracks.csv <- function(file, ...) {
	data.raw <- read.csv(file, ...)
	if (ncol(data.raw) == 4) {
		colnames(data.raw) <- c("id", "pos", "t", "x")
	}
	if (ncol(data.raw) == 5) {
		colnames(data.raw) <- c("id", "pos", "t", "x", "y")
	}
	if (ncol(data.raw) == 6) {
		colnames(data.raw) <- c("id", "pos", "t", "x", "y", "z")
	}
	sort.tracks(structure( 
		split( data.raw[,3:ncol(data.raw)], data.raw$id ),
		class="tracks"))
}



#' Plot in 2D Space
#' 
#' Plots tracks contained in a "tracks" object into a twodimensional space
#' pallelel to the data's axes.
#' 
#' @param x the tracks that shall be plotted, in the form of a object of 
#' class "tracks", as returned by e.g. \code{\link{read.tracks.csv}}.
#' @param dims a vector giving the dimensions of the track data that shall be 
#' plotted, e.g. \code{c('x','y')} for the \eqn{x} and \eqn{y} dimension.
#' @param add boolean value indicating whether the tracks are to be added to the
#' current plot.
#' @param col a specification of the color(s) to be used. This can be a vector
#'  of size \code{length(x)}, where each entry specififes the color for the 
#'  corresponding track.
#' @param xlab a sting that should be used as a label for the x axis.
#' @param ylab a sting that should be used as a label for the y axis.
#' @param ... additional parameters to be passed to \code{\link[graphics]{plot}} 
#' (in case add=False) or \code{\link[graphics]{points}} (in case add=True), 
#' respectively.
#' @details One dimension of the data (by default \eqn{y}) is plotted against 
#' another (by default \eqn{x}). The dimesions can be chosen by means of the 
#' parameter \code{dims} and the axes can be labeled accordingly with the aid 
#' of \code{xlab} and \code{ylab}. The color can be set through \code{col}.
#' If the tracks should be added to an existing plot, \code{add} is to be set 
#' to \code{TRUE}.
plot.tracks <- function(x, dims=c('x','y'), add=F,
		col=order(names(x)), xlab="X position", ylab="Y position", ... ) {
	args <- list(...)

	ids <- names(x)
  	lcol <- rep_len(col, length(ids))
  	pcol <- rep(lcol, lapply(x,nrow))
	px <- .ulapply(x,'[[',dims[1])
	py <- .ulapply(x,'[[',dims[2])

	if (add==T) {
		points( px, py, col=pcol, cex=.5, ... )
	} else {
		plot( px, py, xlab=xlab, ylab=ylab, col=pcol, cex=.5, ...)
	}
  
	starting.points <- lapply( x, function(p) p[1,] )
	px <- sapply(starting.points,'[[',dims[1])
	py <- sapply(starting.points,'[[',dims[2])
	points( px, py, col=col, cex=2 )
  
	lseg <- function(f,i) {
	    function(d) f( d[,dims[i]], -1 ) 
	}

	x0 <- .ulapply(x,lseg(head,1))
	y0 <- .ulapply(x,lseg(head,2))

	x1 <- .ulapply(x,lseg(tail,1))
	y1 <- .ulapply(x,lseg(tail,2))

	lcol <- rep(lcol, lapply(x, function(d) {
      nrow(d) - 1
    }))
	segments(x0, y0, x1, y1, col=lcol)
}

#' Wrap a Single Track into a \emph{Tracks} Object
#' 
#' Makes a \emph{tracks} object out of a single track.
#' 
#' @param track the track that is to be wrapped.
#' 
#' @details A list containing only the input \code{track} is returned. It is 
#' assigned the name "1". This list is returned as a \emph{tracks} 
#' object.
#'
#' @return Retruns a \emph{tracks} object containing the input \code{track}
#' named "1".
wrapTrack <- function(track) {
  return(as.tracks(list("1" = track)))
}

# TODO handle text duplication in details and return 
#' Staggered Version of a Function
#' 
#' Returns the staggered version of a track measure.
#' 
#' @param measure the measure that shall be computed in a staggered fashion.
#' @param ... further parameters for staggered computation of the measure 
#' (see \code{\link{computeStaggered}}).
#' 
#' @details Returns a funtion of a cell track (in a data frame) that computes
#' the given measure in a staggered fashion on that track.
#' 
#' @return Returns a funtion of a cell track (in a data frame) that computes
#' the given measure in a staggered fashion on that track.
staggered <- function(measure, ...){
  return(function(track) { 
    computeStaggered(track, measure, ...) 
  })
}


# TODO make examples (last line of last param) a link?
#' Compute a Measure on a Track in a Staggered Fashion
#' 
#' Computes a measure on all a track's subtracks and either returns them in a 
#' matrix or returns their mean.
#' 
#' @param track the track for which the measure is to be computed.
#' @param measure the measure that is to be computed.
#' @param matrix a logical indicating whether the whole matrix of values for 
#' the measure for each of the input track's subtracks is to be returned. 
#' Otherwise only the mean is returned.
#' @param min.segments the number of segments that each regarded subtrack 
#' should at least consist of. Typically, this value would be set to the 
#' minimum number of segments that a (sub)track must have in order for the 
#' measure to be decently computed. For example, for the funtion 
#' \code{\link{overallAngle}} at least two segments are necessary. The default
#' value is 1.     
#' 
#' @details The measure is computed for each of the input track's subtracks of 
#' length at least \code{min.segments}, and the resulting values are either 
#' returned in a matrix (if \code{matrix} is set), or their mean is returned. 
#' Note that the returned matix is symmetric, because the direction in which a 
#' track is passed is assumed not to matter, and the values at 
#' \code{[i, i + j]}, where j is a nonnegative integer with 
#' \eqn{j < }\code{min.segments}, (with the default value \code{min.segments=1}
#' this is exactly the main diagonal) are set to \code{NA}.
#'  
#' @return If \code{matrix} is set, a matrix with the values of the measure for
#' all the input track's subtracks is returned. The value of this matrix at 
#' position \code{[i, j]} corresponds to the subtrack that starts with the input track's
#' \eqn{j}th point and ends at its \eqn{i}th. Thus, with increasing column number,
#' the regarded subtrack's starting point is advanced on the original track, 
#' and the values for increasingly long subtracks starting from this point can 
#' be found columnwise below the main diagonal, responsively.
#' If \code{matrix} is not set, the mean over the values of the measure for all
#' subtracks of at least \code{min.segments} segments is retruned.
#' 
#' @examples 
#' # Compute the staggered matrix for overallAngle appilied to all adecuate 
#' # subtracks of the first TCell track
#' computeStaggered(getTracks(TCells, "1")$'1', overallAngle, matrix=TRUE, 
#' min.segments = 2)
computeStaggered <- function(track, measure, matrix=FALSE, min.segments=1) {
  ## old method slightly slower (ca. 9 s vs. 10 s in 1000 replications)
  lag.track <- data.frame(track)
  if (matrix) {
		mat <- matrix(nrow=nrow(track), ncol=nrow(track), 0)
    diag(mat) <- NA
  } else {
    stag.meas <- c()
  }
  if (!matrix) {
    segMeans <- .computeSegmentwiseMeans(track, measure, min.segments)  
    n.before <- length(segMeans)
    segMeans <- segMeans[!is.nan(segMeans)]   # es entstehen NaNs, wenn die Zelle mehrere Zeitschritte an derselben stelle bleibt -> diese werden ignoriert; Parameter?
    n <- length(segMeans)
    if (n.before!=n) {
      warning("NaNs have been removed.")
    }
    weights <- seq(n, 1)
    weights <- weights[!is.nan(segMeans)]
    ret <- weighted.mean(segMeans, weights)
    return(ret)
  } else {
    for (i in (min.segments):(nrow(track) - 2)) {
      subtracks <- computeSubtracksOfISegments(track, i)
      val <- sapply(subtracks, measure)
      diag(mat[-1:-i,]) <- val
    }
    mat[nrow(track), 1] <- sapply(computeSubtracksOfISegments(track, (nrow(track)-1)), measure)
    mat <- mat + t(mat)
    return(mat)
    }                                        
}


# TODO handle text duplication in details and return
#' Measure for Every Prefix of a Track
#' 
#' Returns a function that applies the input measure to each of a tracks prefixes 
#' of a given minimal length.
#' 
#' @param measure the measure that is to be computed.
#' @param min.length minimal number of points in a prefix the measure is to be 
#' applied to. Default is \eqn{1}.
#' 
#' @details A function that applies the measure to every prefix of length at 
#' least \code{min.length} is returned. This 
#' function expects a track as an input, whose prefixes are considered.
#' 
#' @return A function that applies the measure to every prefix of length at 
#' least \code{min.length} is returned. This 
#' function expects a track as an input, whose prefixes are considered.
forEveryPrefix <- function(measure, min.length=1) {
  function(track) {
    computeForEveryPrefix(track, measure, min.length) 
  }
}

# formerName: computeForEveryTruePrefixOfLengthAtLeast
#' Apply a Measure to Every Prefix of a Given Track
#' 
#' Computes a vector in which each position represents the value of the given 
#' measure for a prefix of certain length of the input track.
#' 
#' @param track the track whose prefixes are to be regarded.
#' @param measure the measure that is to be applied.
#' @param min.length minimal number of points in regarded prefixes. Default is 
#' \eqn{1}.
#' 
#' @details The measure is applied to every prefix of the input track, starting
#' with the prefix that consists of \code{min.length} points. The resulting 
#' values are returned in a vector, sorted in ascending order of the 
#' corresponding prefix's length.
#' 
#' @return Returns a vector with the values of the measure applied to each of 
#' the given track's prefixes.
computeForEveryPrefix <- function(track, measure, min.length=1) {
  res <- c()
  for (i in (min.length:nrow(track))) {
    res <- c(res, measure(track[1:i,]))
  }
  return(res)
}


#' Normalize a track
#' 
#' Translates a track such that its starting point is in the origin of ordinates.
#' 
#' @param track the track that is to be translated, given in a data frame with 
#' columns \eqn{t}, \eqn{x} and possibly \eqn{y}, and \eqn{z}.
#' 
#' @details Subtracts the coodinates and time of the tracks starting point from
#' every point on the track.
#' 
#' @return Retruns the translation of the input track whose starting point is 
#' on the origin of ordinates.
normalizeTrack <- function(track) {
  # Subtract first row from each row
  norm.track <- data.frame(track$t, 
                           Reduce(rbind, apply(as.matrix(track[, 2:ncol(track)]), 
                                               1, function(v) {
                                                 v - track[1, 2:ncol(track)]
                                               }), c()))
  colnames(norm.track) <- colnames(track)
  return(norm.track)
}


#' Normalize Tracks
#' 
#' Normalizes list of tracks (given in a \code{tracks} object) such that the 
#' starting point of each track is in the origin of ordinates.
#' 
#' @param tracks the \code{tracks} object that is to be normalized.
#' 
#' @details The function \code{\link{normalizeTrack}} is applied to each track.
#' 
#' @return A \emph{tracks} object with the tracks from the input object translated 
#' such that their starting point is in the origin of ordinates.
normalizeTracks <- function(tracks) as.tracks(lapply(tracks, normalizeTrack))




#' All Subtracks of a Given Length
#' 
#' For a given integer \eqn{i} and a track, this function returns a tracks 
#' object of all subtracks that consist of exactly \eqn{i} subsegments. 
#' 
#' @param track the track whose subtracks shall be computed.
#' @param i integer indicating the desired subtrack length
#' @param overlap the number of segments that shall overlap in each two 
#' subsequent subtracks.
#' 
#' @details Beginning with the track's starting point, for each point on
#' the track, this point and the subsequent \eqn{i} points are saved 
#' in a data frame. All the data frames are returned as a list, constituting
#' an object of class \emph{tracks}.
#' 
#' @return A list of data frames (a \emph{tracks} object), containing a data
#' frame for each of the track's subtracks that constists of exactly \eqn{i} single 
#' segments.
computeSubtracksOfISegments <- function(track, i, overlap=i-1) {
  if (nrow(track) <= i) {
    return(NA)
  } 
  if (overlap > i-1) {
    warning("Overlap exceeds segment length")
    overlap <- i-1
  }
  l <- list()
#   print(paste("i:", i, "overlap", overlap))
  for (j in seq(1, (nrow(track)-i),i-overlap)) {
#     print(paste(j, j+i))
    l[[as.character(j)]] <- track[j:(j+i), ]    #l[[as.character(j)]]
  }    
  structure(l, class="tracks")
}


#' All Subtracks with a Certain Number of Segments for a \emph{Tracks} Object
#' 
#' Computes a \emph{tracks} object of all the subtracks of the input 
#' \emph{tracks} object that constist of \code{i} segments.
#' 
#' @param tracks the \emph{tracks} object whose subtracks are to be computed.
#' @param i the number of segments the subtracks ahll consist of.
#' @param overlap the number of segments in which each subtrack shall overlap 
#' its predecessor and successor. Default is \eqn{i - 1}.
#' 
#' @details The function \code{\link{computeSubtracksOfISegments}} is applied 
#' on every track for the given \code{i} and \code{overlap} and the result is 
#' returned as a single \emph{tracks} object. 
#' 
#' @return A \emph{tracks} object is returned which contains all the subtracks 
#' of any track in the input \emph{tracks} object that consist of exactly \code{i}
#' segments and overlap their predecessor in  \code{overlap} segments.
subtracks <- function(tracks, i, overlap=i-1 ) {
  Reduce(c, lapply(tracks, subtracksOfISegments(i, overlap)))
}

# TODO remove test duplications in descr., details and return, merge information with above function
#' Measure for Every Subtrack of a Track
#' 
#' Returns a function that for a given track outputs all subtracks that 
#' constist of a certain number of segments. The overlap of subsequent 
#' subtracks can also be given.
#' 
#' @param i the number of segments the outputted subtracks shall consist of.
#' @param overlap number no segments that shall be common to subsequent 
#' subtracks. Default is \eqn{i-1}.
#' 
#' @details Returns a function that for a given track outputs all subtracks that 
#' constist of a certain number of segments. The outputted segments are sorted 
#' by their starting point (the one starting with the input track's starting 
#' point will be the first), and two subsequent subtracks will have 
#' \code{overlap} segments in common, by default \eqn{i-1}, so they differ in 
#' only one segment.
#'  
#' @return Returns a function that for a given track outputs all subtracks that 
#' constist of a certain number of segments. The outputted segments are sorted 
#' by their starting point (the one starting with the input track's starting 
#' point will be the first), and two subsequent subtracks will have 
#' \code{overlap} segments in common, by default \eqn{i-1}, so they differ in 
#' only one segment.
subtracksOfISegments <- function(i, overlap=i-1) {
  function(track) {
    computeSubtracksOfISegments(track, i, overlap)
  }
}


#' Normalize a Measure to Track Duration
#' 
#' Returns a measure that divides the input measure by the duration of its 
#' input track.
#' 
#' @param measure the measure that should be nomalized.
#' 
#' @details A function is retruned that computes the input measure for a given 
#' track an divides the result by the track's duration.
#' 
#' @return Retruns a function that computes the input measure for a given track
#' and returns the result divided by the track's duration.
#' 
#' @examples # normalizeToDuration(displacement) can be used as an indicator 
#' # for the motions efficiency
#' lapply(TCells, normalizeToDuration(displacement))
normalizeToDuration <- function(measure) {
  function(track) {
    measure(track) / duration(track)
  }
}



#' Bivariate Scatterplot of 2 Given Track Measures
#' 
#' Plots the values of two measures applied on given tracks against each 
#' other.
#' 
#' @param measurex the measure to be shown on the X axis.
#' @param measurey the measure to be shown on the Y axis.
#' @param x the argument which the measures shall be applied to. Either this 
#' must be of a type which \code{measurex} and \code{measurey} can directly
#' be applied to, or, if measurex and measurey are applicable to single tracks,
#' the argument can be a \emph{tracks} object, which the measures will be 
#' applied to using sapply before plotting the result.
#' @param add a logical indicating whether the tracks are to be added to an 
#' existing plot via \code{\link[graphics]{points}}.
#' @param xlab label of the x-axis. By default the name of the input function 
#' \code{measurex}.
#' @param ylab label of the y-axis. By default the name of the input function 
#' \code{measurey}.
#' @param ... additional parameters to be passed to \code{\link[graphics]{plot}} 
#' (in case add=False) or \code{\link[graphics]{points}} (in case add=True), 
#' respectively.
#' 
#' @details Plots the value of \code{measurey} applied to \code{x} against the
#' value of  \code{measurey} applied to \code{y}.
plotTrackMeasures <- function(measurex, measurey, x, add=FALSE, 
                              xlab=deparse(substitute(measurex)), 
                              ylab=deparse(substitute(measurey)), ...) {
  if(class(x)=="tracks") {
    if(add) {
      points(sapply(x,measurex) ~ sapply(x,measurey), ...)
    } else {
      plot(sapply(x,measurex) ~ sapply(x,measurey), xlab=xlab, ylab=ylab, ...)
    }
  } else {
    if(add) {
      points(measurex(x) ~ measurey(x), ...)
    } else {
      plot(measurex(x) ~ measurey(x), xlab=xlab, ylab=ylab, ...)
    }
  }
}
  

#' Cluster Tracks
#' 
#' Clusters the tracks in a given \emph{tracks} object according to given 
#' measures.
#' 
#' @param tracks the tracks that are to be clustered.
#' @param measures the measures according to which the tracks are to be 
#' clustered.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the clustering. Default 
#' is TRUE.
#' @param ... additional parameters to be passed to \code{\link[stats]{hclust}}.
#' 
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object and according to the resulting values the tracks are 
#' clustered using a hierarchical clustering (see \code{\link[stats]{hclust}}).
#' If \code{scale} is set to true, the values that have been determined by 
#' applying the input measures to the input tracks are scaled to mean value 
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the clustering.
clusterTracks <- function(tracks, measures, scale=TRUE, ...) { 
	values <- matrix(nrow=length(tracks))
	if (is.function(measures)) {
		measures <- c(measures)
	}
	values <- Reduce(cbind, lapply(measures, function(m) sapply(tracks, m)))
	if (scale) {
		values <- scale(values)
	}
	hclust(dist(values), ...)
}


# TODO: Warum wird die Null im defaultwert fuer measurey so nach hinten verschoben? (schon im Rd-file)
#list version bringt wegen skalierung der Achsen nichts 
#' Plot Tracks in a 1D or 2D Parameter Space
#' 
#' Plots the tracks in a \emph{tracks} object into the parameter space of one 
#' or two given measures.
#' 
#' @param tracks a \emph{tracks} object containing the tracks that are to be 
#' plotted
#' @param measurex the measure which constitutes the x-dimension of the 
#' parameter space
#' @param measurey the measure which constitutes the y-dimension of the 
#' parameter space. By default constant \eqn{0}.
#' @param add a logical indicating whether the tracks are to be added to an 
#' existing plot via \code{\link[graphics]{points}}.
#' @param xlab label of the x-axis. By default the name of the input function 
#' \code{measurex}.
#' @param ylab label of the y-axis. By default the name of the input function 
#' \code{measurey}.
#' @param col the color(s) that are to be used, by default color 1 (black), 
#' see \code{\link[graphics]{par}}
#' @param ... additional parameters to be passed to \code{\link[graphics]{plot}} 
#' (in case add=False) or \code{\link[graphics]{points}} (in case add=True), 
#' respectively.
#' 
#' 
#' @details Plots the tracks in the given \emph{tracks} object into the 
#' parameter space of \code{measurex} and \code{measurey}, or of 
#' \code{measurex} and \eqn{0}, if \code{measurey} is not given. 
#' 
plotInParameterSpace <- function(tracks, measurex, measurey=function(x) {0}, 
                                 add=FALSE, xlab=deparse(substitute(measurex)),
                                 ylab=deparse(substitute(measurey)), col=1, ...) {
    if (add) {
      points(unlist(lapply(tracks, measurey)) ~ unlist(lapply(tracks, measurex)),  col=col, ...)
    } else {
      plot(unlist(lapply(tracks, measurey)) ~ unlist(lapply(tracks, measurex)), 
           xlab=xlab, ylab=ylab, col=col, ...)
    }
}



  
#' Compute a Statistic on Measure Values of Subtracks of a Certain Length
#' 
#' For a given \emph{tracks} object, a given measure is applied to all the 
#' tracks' subtracks of a certain number of segments \eqn{i}. 
#' On all the values obtained for a certain value of \eqn{i}, the given 
#' statistic is computed, and all the resulting values are returned.
#' 
#' @param tracks the tracks object whose subtracks are to be regarded.
#' @param number.of.segments either a numeric value \eqn{i}, indicating that only 
#' subtracks of exactly \eqn{i} segments are to be regarded, or a vector of all
#' the values \eqn{i} that shall be used.
#' @param measure the measure that is to be computed on the subtracks.
#' @param statistic the function that is to be used to aggregate the measures 
#' of subtracks of the same number of segments. The function may return an 
#' arbitrary number of values.
#' For the most common functions, a string can be given:
#' \itemize{
#'  \item{"mean.sd"}{ Outputs the mean and \eqn{mean - sd} as lower and 
#'  \eqn{mean + sd} as upper bound}
#'  \item{"mean.se"}{ Outputs the mean and \eqn{mean - se} as lower and 
#'  \eqn{mean + se} as upper bound}
#'  \item{"mean.ci.95"}{ Outputs the mean and upper and lower bound of the 
#'  symmetric 95 percent confidence intervall}. 
#'  \item{"mean.ci.99"}{ Outputs the mean and upper and lower bound of the 
#'  symmetric 95 percent confidence intervall}. 
#'  \item{"iqr"}{ Outputs the interquartil range, that is, the median, and the 
#'  25-percent-quartile as a lower and and the 75-percent-quartile as an 
#'  upper bound}
#' }
#' @details For every number of segments \eqn{i} in the set defined by 
#' \code{number.of.segments}, all subtracks of any track in the input 
#' \emph{\code{tracks}} object, that consist of exactly \eqn{i} segments are 
#' regarded. The input \code{measure} is applied to the subtracks individually, 
#' and the \code{statistic} is used to agglomerate the obtained values. 
#' The return values of the agglomeration function are returned for each value 
#' of \eqn{i} in a data frame.
#' @return Returns a data frame with one line for every \eqn{i} in the set 
#' specified by \code{number.of.segments}. The first column contains the values 
#' of \eqn{i}, the following columns contain the values of the agglomeration 
#' function applied on the measure values of tracks of exactly \eqn{i} segments.
#' @examples 
#' require(ggplot2)
#' dat <- aggregateSubtrackMeasures(TCells, displacement, "mean.se")
#' ggplot(dat, aes(x=i, y=mean)) + geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)
aggregateSubtrackMeasures <- function(tracks, measure, statistic=mean, 
                                      number.of.segments=seq(1, (maxTrackLength(tracks)-1))) {  # geht fuer einzelne tracks nicht!
#   if(class(number.of.segments)=="numeric") {
#     subtracks <- computeSubtracksOfISegments(track, number.of.segments) ## overlap uebergeben???
#   } else {
  if (max(number.of.segments) > (maxTrackLength(tracks)-1)) {
    warning("No track is long enough!")
  }
  if (is.character(statistic)) {
    if (statistic == "mean.ci.95") {
      statistic <- function(x) {
        ci <- tryCatch( t.test(x)$conf.int, error=function(e) rep(mean(x),2) )
        return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
      }
    } else if (statistic == "mean.ci.99") {
      statistic <- function(x) {
        ci <- tryCatch( t.test(x, conf.level=.99)$conf.int, error=function(e) rep(mean(x),2) )
        return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
      }
    } else if (statistic == "mean.se") {
      statistic <- function(x) {
        return(c(mean=mean(x), lower = mean(x) - sd(x)/sqrt(length(x)), 
                 upper = mean(x) + sd(x)/sqrt(length(x))))
      }
    } else if (statistic == "mean.sd") {
      statistic <- function(x) {
        return(c(mean=mean(x), lower = mean(x) - sd(x), 
                 upper = mean(x) + sd(x)))
      }
    } else if (statistic == "iqr") {
      statistic <- function(x) {
        lx <- quantile(x,probs=c(.25,.5,.75))
        names(lx) <- c("lower","median","upper")
        return(lx)
      }
    } else {
      stop("Statistic not known")
    }
  } else {
    if (!is.function(statistic)) {
       stop("statistic must be a function or a string")
    }
  }
  subtracks <- list()
  measure.values <- list()
  for (i in number.of.segments) {
    subtracks[[i]] <- subtracks(tracks, i)   ## overlap uebergeben???
    subtracks[[i]] <- subtracks[[i]][sapply(subtracks[[i]], function(z) !all(is.na(z)))]     # remove NULL elements from list
    if(class(subtracks[[i]]) == "list" ) {
      measure.values[[i]] <- sapply(subtracks[[i]], measure)
    } else {
      measure.values[[i]] <- measure(subtracks[[i]])
    }
  }
  measure.values <- measure.values[!sapply(measure.values, is.null)] # remove NULL elements from list
  value <- sapply(measure.values, statistic)
#   names(ret) <- deparse(substitute(statistic))  ##wie bekommt man bei einer eindimensionalen Fkt. die Spalte nach dem Fkt.namen benannt?
#   stat = deparse(substitute(statistic))
  ret <- rbind(i=number.of.segments, value)
  ret <- t(ret)
  return(data.frame(ret))
}



#' The Maximum Number of Points of a Track in a \emph{Tracks} Object
#' 
#' Determines the maximum number of points that one of the tracks in the given
#' \emph{tracks} object constists of.
#' 
#' @param tracks the \emph{tracks} object the tracks in which are to be considered.
#' 
#' @details The maximum number of rows of a track in the \emph{tracks} object 
#' is determined.
#' 
#' @return Returns the maximum number of rows of a track in the given 
#' \emph{tracks} object.
maxTrackLength <- function(tracks) {
  return(max(sapply(tracks, nrow)))
}



# TODO remove for-loop
#' The Bounding Box of a \emph{Tracks} Object
#' 
#' Computes minimum and maximum value that a track has in each dimension
#' 
#' @param tracks the \emph{tracks} object whose bounding box is to be computed
#' 
#' @details In each of the tracks' dimensions the minimum and maximum value 
#' that any track from the input has is computed and returned. It is assumed 
#' that all the tracks have the same number of dimensions. 
#'  
#' @return Returns a matrix with two rows and \eqn{d} columns, where \eqn{d} is 
#' the number of dimensions of the tracks. The first row contains the minimum 
#' and the second row the maximum value of any track in the dimension given by 
#' the column.
boundingBox <- function(tracks) {
  if( class(tracks) != "tracks" ){
    tracks <- as.tracks( list(T1=tracks) ) 
  }
  bounding.boxes <- lapply(tracks, .boundingBoxTrack)
#   print(bounding.boxes)
 
  empty <- mat.or.vec(nrow(bounding.boxes[[1]]), ncol(bounding.boxes[[1]]))
  colnames(empty) <- colnames(bounding.boxes[[1]])
  rownames(empty) <- rownames(bounding.boxes[[1]])
  empty[1,] <- Inf
  empty[2,] <- -Inf
#   print(empty)
  return(Reduce(.minMaxOfTwoMatrices, bounding.boxes, empty)) 
}

# TODO ellipse einzeichnen?
#' Unbiasedness of Motion
#' 
#' Determines the p-value for the motion described by a \emph{tracks} object to 
#' be unbiased.
#' 
#' @param tracks the tracks whose biasedness is to be determined.
#' @param dim vector with the names of the track's 
#' dimensions that are to be considered. By default c("x", "y").
#' @param plot logical indicating whether the scatter of the step's directions, 
#' origin of ordinates (green circle) and the mean of the data points (green 
#' cross) are to be plotted. (In one dimension also the bounds of the 
#' condfidence interval are given.) Plot works only in one or two dimensions.
#' 
#' @details Computes the displacement vectors of all segments in the tracks 
#' given in \code{tracks}, and determines the p-value for this distribution 
#' using Hotelling's T-square Test 
#' (see \code{\link[DescTools]{HotellingsT2Test}}).
#' 
#' 
hotellingsTest <- function(tracks, dim=c("x", "y"), plot=FALSE) {
  stopifnot( !plot || (length(dim) %in% c(1,2) ))
  Tx <- projectDimensions(tracks, dim)
  if (length(dim) > 1) {
    sx <- t(Reduce(cbind,(lapply(1:length(Tx), 
                                 function(i) sapply(subtracks(Tx[i], 1),
                                                    function(x) unlist(displacementVector(x)))))))
    
    if (plot) {
      plot(sx)
      points(0, 0, col=3)
      points(mean(sx[, 1:ncol(sx)]), col=2, pch=4, cex=2)
    }
  } else {
    sx <- Reduce(c,(lapply(1:length(Tx), 
                                 function(i) sapply(subtracks(Tx[i], 1),
                                                    function(x) unlist(displacementVector(x))))))
    if (plot) {
      mean.sx <- mean(sx)
      sd.sx <- sd(sx)
      plot(0, 0, col=3, xlab = dim[1])
      stripchart(sx, add=TRUE, at=0, pch=1)
      points(mean.sx, 0, col=2, pch=4, cex=2)
      points(0, 0, col=3)
      print(mean.sx -  sd.sx)
      points(mean.sx - sd.sx / sqrt(length(sx)), 0, pch=3, col=4)
      points(mean.sx + sd.sx / sqrt(length(sx)), 0, pch=3, col=4)
    }
  }
  return(DescTools::HotellingsT2Test(sx)[["p.value"]])
}

byPrefixLength <- function( tracks, f, simplify=TRUE ){
  
}