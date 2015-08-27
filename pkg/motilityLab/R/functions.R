#' Convert Tracks to Data Frame
#' 
#' Converts a \code{tracks} object into a data frame.
#' 
#' @param x the \code{tracks} object to be coerced to a data frame.
#' 
#' @param row.names NULL or a character vector giving row names for the
#'  data frame.  Missing values are not allowed.
#' @param optional logical. Required for S3 consistency, but 
#' has no effect: column names are always assigned to the resulting
#'  data frame regardless of the setting of this option.
#' @param ... further arguments to be passed from or to other methods.
#'
#' @details Converts tracks from the list-of-matrices format, which is good
#' for efficient processing and therefore the default in this package, to a 
#' single dataframe which is efficient for tasks such as plotting or 
#' writing a dataset to a CSV file. 
#"
#' @return a single data frame containing all individual tracks from the input with a 
#' prepended column named "id" containing each track's identifier in `x`.
as.data.frame.tracks <- function(x, row.names = NULL, optional = FALSE, ...) {
	ids <- rep(names(x), lapply(x,nrow))
	r <- data.frame(id=ids,  Reduce( function(...) 
		rbind.data.frame(...,make.row.names=FALSE), x ) )
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
#' @param ... further arguments to be passed on to the respective method.
#'
#' @details The S3 Method corresponding to the class of `x` (usually \code{list})
#' is called.
#'  
#' @return the coerced object of S3 class \code{tracks}.
as.tracks <- function(x,...) 
  UseMethod("as.tracks")


#' Convert from Tracks to List
#' 
#' Coerces a \code{tracks} object to a list
#' 
#' @param x the \code{tracks} object to be coerced to a list.
#' @param ... further arguments to be passed from or to other methods.
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
#' @param ... further arguments to be passed from or to other methods.
#' 
#' @details The input list is assigned the S3 class \code{tracks}.
#' 
#' @return the coerced object of S3 class \code{tracks} (i.e., a list 
#'  containing a data frame for each track)
as.tracks.list <- function(x,...) 
  structure(x, class="tracks")


#' Get a Subset of Tracks
#' 
#' Returns a \code{tracks} object containing only the specified elements of a given
#' \code{tracks} object.
#' 
#' @param x the input tracks object.
#' @param ids a vector containing the ids of the tracks that shall be chosen.
#' 
#' @details the elements given in \code{ids} are chosen from the list that represents
#' the \code{tracks} object and returned, likewise as a \code{tracks} object.
#' 
#' @return the modified \code{tracks} object.
getTracks <- function(x, ids) {
	structure(x[ids], class="tracks")
}

#' Remove Tracks from Tracks Object
#' 
#' Returns a \emph{tracks} object from the tracks with the given IDs are removed.
#'
#' @param x the input tracks object.
#' @param ids a vector containing the ids of the tracks that shall be removed.
#'
#' @return the modified \code{tracks} object.
removeTracks <- function(x,ids) {
	structure(x[!(names(x) %in% ids)], class="tracks")
}

#' Filter Tracks 
#'
#' Extracts subtracks based on a given function.
#'
#' @param f a function that accepts a single track as its first argument and returns a 
#' logical value (or a value that can be coerced to a locical).
#' @param x a tracks object.
#' @param ... further arguments to be passed on to \code{f}.
#' @return a \code{tracks} object containing only those tracks from \code{x} for which
#' \code{f} evaluates to \code{TRUE}.
#' @examples
#' ## Remove short tracks from the T cells data
#' plot( filterTracks( function(t) nrow(t)>10, TCells ) )
filterTracks <- function(f,x,...){
	structure(x[as.logical(sapply(x,function(x) f(x,...) ))], class="tracks")
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
#' @param ... further arguments to be passed on to \code{order}.
#'
#' @details Sorts the positions of each track (represented as a data frame) in the 
#' \emph{tracks} object by time (given in the column t). 
#' 
#' @return A \emph{tracks object} that contains the tracks from the input object 
#' sorted by time is returned.
sort.tracks <- function(x, decreasing=FALSE, ...) {
	as.tracks(lapply(x, function(d) d[order(d[,"t"],decreasing=decreasing),] ))
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


#' Read Tracks from CSV File
#' 
#' Reads cell tracks from a CSV file. Data are expected to be organized as follows.
#' One column contains a track identifier, which can be numeric or a string, and 
#' determines which points belong to the same track. 
#' Another column is expected to contain a time index or a time period (e.g. number of
#' seconds elapsed since the beginning of the track, or since the beginning of the 
#' experiment). Input of dates is not (yet) supported, as absolute time information is
#' frequently not available. 
#' One to three further columns contain the spatial coordinates 
#' (depending on whether the tracks are 1D, 2D or 3D). 
#' The names or indices of these columns in the CSV files are given using the 
#' corresponding parameters (see below). Names and indices can be mixed, e.g. you can
#' specify \code{id.column="Parent"} and \code{pos.columns=1:3}
#' 
#' @param file the name of the file which the data are to be read from, a 
#' readable text-mode connection or a complete URL
#' (see \code{\link[utils]{read.table}}).
#'
#' @param id.column index or name of the column that contains the track ID.
#'
#' @param time.column index or name of the column that contains elapsed time.
#'
#' @param pos.columns vector containing indices or names of the columns that contain
#' the spatial coordinates.
#'
#' @param scale.t a value by which to multiply each time point. Useful for changing units,
#' or for specifying the time between positions if this is not contained in the file
#' itself.
#'
#' @param scale.pos a value, or a vector of values, by which to multiply each spatial 
#' position. Useful for changing units.
#' 
#' @param ... Further arguments to be passed to \code{read.csv}, for instance 
#' \code{sep="\t"} can be useful for tab-separated files.
#' 
#' @return An object of class \emph{tracks} is returned, which is a list of 
#' matrices, each containing the positions of one track. The matrices 
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
read.tracks.csv <- function(file, id.column=1, time.column=2, 
	pos.columns=c(3,4,5), scale.t=1, scale.pos=1, ...) {
	data.raw <- read.csv(file, ...)
	if( ncol(data.raw) < length(pos.columns)+2 ){
		stop("CSV file does not contain enough columns! (Perhaps you need to specify 'sep')")
	}
	if( length(pos.columns) < 1 ){
		stop("At least one position column needs to be specified!")
	}
	cx <- as.character(c(id.column,time.column,pos.columns))
	cxc <- match( cx, colnames(data.raw) )
	cxi <- match( cx, seq_len(ncol(data.raw)) )

	cxi[is.na(cxi)] <- cxc[is.na(cxi)]

	if( any(is.na(cxi)) ){
		stop("Column(s) not found: ",
			paste(cx[is.na(cxi)],collapse=","))
	}
	r <- data.raw[,as.integer(cxi)]
	if( ncol(r) <= 5 ){
		colnames(r) <- c("id","t",c("x","y","z")[seq_along(pos.columns)])
	} else {
		colnames(r) <- c("id","t",paste0("x",seq_along(pos.columns)))
	}
	if( scale.t != 1 ){
		r[,"t"] <- scale.t*r[,"t"]
	}
	if( any( scale.pos != 1 ) ){
		r[,-c(1,2)] <- scale.t*r[,-c(1,2)]
	}
	sort.tracks(structure( 
		split.data.frame( as.matrix(r[,2:ncol(r)]), r[,1] ),
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
	px <- unlist( lapply(x,'[', , dims[1]) )
	py <- unlist( lapply(x,'[', , dims[2]) )

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

#' Staggered Version of a Function
#' 
#' Returns the "staggered" version of a track measure. That is, instead of
#' computing the measure on the whole track, the measure is average over
#' all subtracks (of any length) of the track.
#' 
#' @param measure the measure that shall be computed in a staggered fashion.
#' @param ... further parameters passed on to \code{\link{computeStaggered}}.
#' 
#' @details This is a wrapper mainly designed to provide a convenient interface
#' for track-based staggered computations with \code{lapply}, see example.
#' 
#' @return Returns a function that computes
#' the given measure in a staggered fashion on that track.
#' 
#' @references \cite{Mokhtari et al.: Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior, PLoS ONE 2013 Dec 6;8(12):e80808, 2013} 
#'
#' @examples
#' hist( sapply( TCells, staggered( displacement ) ) )
#'
staggered <- function(measure, ...){
  return( function(track) computeStaggered(track, measure, ...) )
}


#' Compute a Measure on a Track in a Staggered Fashion
#' 
#' Computes a measure on all subtracks of a track and return them either
#' as a matrix or return their mean.
#' 
#' @param track the track for which the measure is to be computed.
#' @param measure the measure that is to be computed.
#' @param matrix a logical indicating whether the whole matrix of values for 
#' the measure for each of the input track's subtracks is to be returned. 
#' Otherwise only the mean is returned.
#' @param min.segments the number of segments that each regarded subtrack 
#' should at least consist of. Typically, this value would be set to the 
#' minimum number of segments that a (sub)track must have in order for the 
#' measure to be decently computed. For example, at least two segments are needed
#' to compute the \code{\link{overallAngle}}.
#' 
#' @details The measure is computed for each of the input track's subtracks of 
#' length at least \code{min.segments}, and the resulting values are either 
#' returned in a matrix (if \code{matrix} is set), or their mean is returned. 
#' The computed matrix is symmetric since the direction along which a 
#' track is traversed is assumed not to matter. The values at 
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
#' be found columnwise below the main diagonal, respectively.
#' If `matrix` is not set, the mean over the values of the measure for all
#' subtracks of at least `min.segments` segments is retruned.
#' 
#' @examples 
#' # Compute the staggered matrix for overallAngle applied to all adequate 
#' # subtracks of the first TCell track
#' computeStaggered(TCells[[1]], overallAngle, matrix=TRUE, min.segments = 2)
computeStaggered <- function(track, measure, matrix=FALSE, min.segments=1) {
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
      subtracks <- subtracks(track, i)
      val <- sapply(subtracks, measure)
      diag(mat[-1:-i,]) <- val
    }
    mat[nrow(track), 1] <- sapply(subtracks(track, (nrow(track)-1)), measure)
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
#' function applies the measure to every prefix of the input track, starting
#' with the prefix that consists of \code{min.length} points. The resulting 
#' values are returned in a vector, sorted in ascending order by the 
#' prefix length.
#'
#' @return A function that applies the measure to every prefix of length at 
#' least \code{min.length} is returned. 
#'
forEveryPrefix <- function(measure, min.length=1) {
  function(track) {
	  res <- c()
	  for (i in (min.length:nrow(track))) {
		res <- c(res, measure(track[1:i,,drop=FALSE]))
	  }
	  return(res)
  }
}

#' Normalize a Track
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
	cbind( track[,"t",drop=FALSE], sweep(track[,-1],2,track[1,-1]) )
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

#' Get Subtracks of One or More Tracks
#' 
#' Creates a \emph{tracks} object consisting of all subtracks of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \emph{tracks} object whose 
#' subtracks shall be computed.
#' @param i subtrack length. A single integer, lists are not supported.
#' @param overlap the number of segments in which each subtrack shall overlap 
#' with the previous and next subtrack. The default \code{i - 1} returns all
#' subtracks. Can be negative, which means that space will be left between
#' subtracks. 
#' 
#' @details The output is always a single \emph{tracks} object, which is 
#' convenient for many common analyses. If subtracks are to be considered separately
#' for each track, use the function \code{\link{staggered}} together with
#' \code{lapply}. Subtrack extraction always starts at the first position of the
#' input track.
#' 
#' @return A \emph{tracks} object is returned which contains all the subtracks 
#' of any track in the input \emph{tracks} object that consist of exactly `i`
#' segments and overlap adjacent subtracks in `overlap` segments.
subtracks <- function(x, i, overlap=i-1 ) {
	if( class(x) != "tracks" ){
		x <- wrapTrack( x )
	}
	Reduce(c, lapply(x, 
  		function(t) .computeSubtracksOfISegments(t, i, overlap )))
}

#' Get Track Prefixes 
#' 
#' Creates a \emph{tracks} object consisting of all prefixes (i.e., subtracks
#' starting with the first position of a track) of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \emph{tracks} object whose 
#' prefixes shall be computed.
#' @param i subtrack length. A single integer, lists are not supported.
#'
#' @details This function behaves exactly like \code{\link{subtracks}} except
#' that only subtracks starting from the first position are considered.
#'
prefixes <- function(x,i) {
	if( class(x) != "tracks" ){
		x <- wrapTrack( x )
	}
	lapply( x,
		function(t) {
			if( nrow(t) <= i ){ return( NULL ) }
			t[1:(i+1), ,drop=FALSE]
		}
	)
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
#' @return An object of class *hclust*, see \code{\link[stats]{hclust}}.
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


#' Compute summary statistic of measure values over subtracks of a given track
#' set.
#' 
#' This is the main workhorse function to compute many common motility measures 
#' such as mean square displacement, turning angle, and autocorrelation functions.
#' 
#' For the given \emph{tracks} object, a track measure is applied to all the 
#' tracks' subtracks of a certain number of segments \eqn{i}. 
#' On all the values obtained for a certain value of \eqn{i}, the given 
#' statistic is computed, and the resulting values are returned.
#' 
#' @param x the tracks object whose subtracks are to be considered. 
#' If a single track is given, it will be coerced to a tracks object 
#' using \code{\link{wrapTrack}} (but note that this requires an explicit call
#' \code{aggregate.tracks}).
#' @param measure the measure that is to be computed on the subtracks.
#' @param by a string that indicates how grouping is performed. Currently, two
#' kinds of grouping are supported: 
#' \itemize{
#'  \item{"subtracks"}{ Apply `measure` to all subtracks according to 
#' the parameters `subtrack.length` and `max.overlap`.}
#'  \item{"prefixes"}{ Apply `measure` to all prefixes (i.e., subtracks starting
#'  from a track's initial position) according to the parameter `subtrack.length`.}
#' }
#' 
#' @param FUN a function that is to be used to aggregate the measures 
#' of corresponding subtracks. 
#' If a string is given, it is first matched to the following builtin values:
#' \itemize{
#'  \item{"mean.sd"}{ Outputs the mean and \eqn{mean - sd} as lower and 
#'  \eqn{mean + sd} as upper bound}
#'  \item{"mean.se"}{ Outputs the mean and \eqn{mean - se} as lower and 
#'  \eqn{mean + se} as upper bound}
#'  \item{"mean.ci.95"}{ Outputs the mean and upper and lower bound of the 
#'  symmetric 95 percent confidence intervall}. 
#'  \item{"mean.ci.99"}{ Outputs the mean and upper and lower bound of the 
#'  symmetric 95 percent confidence intervall}. 
#'  \item{"iqr"}{ Outputs the interquartile range, that is, the median, and the 
#'  25-percent-quartile as a lower and and the 75-percent-quartile as an 
#'  upper bound}
#' }
#' If the string is not equal to any of these, it is passed on to 
#' \code{\link{match.fun}}
#'
#' @param subtrack.length either a numeric value \eqn{i}, indicating that only 
#' subtracks of exactly \eqn{i} segments are to be regarded, or a vector of all
#' the values \eqn{i} that shall be used. In particular, subtrack.length=2 corresponds
#' to a "step-based analysis".
#'
#' @param max.overlap Determines the amount by which the subtracks that are taken into
#' account can overlap. A maximum overlap of \code{max(subtrack.length)} will imply
#' that all subtracks are considered. For a maximum overlap of 0, only non-overlapping
#' subtracks are considered. A negative overlap can be used to ensure that only subtracks
#' a certain distance apart are considered. In general, for non-Brownian motion there will
#' be correlations between subsequent steps, such that a negative overlap may be necessary
#' to get a proper error estimate.
#' @param na.rm This is passed on to the builtin statistics like "mean.se" or 
#' "mean".
#' @param ... further arguments passed to or used by methods.
#'
#' @details For every number of segments \eqn{i} in the set defined by 
#' \code{subtrack.length}, all subtracks of any track in the input 
#' \emph{\code{tracks}} object, that consist of exactly \eqn{i} segments are 
#' regarded. The input \code{measure} is applied to the subtracks individually, 
#' and the \code{statistic} is used to agglomerate the obtained values. 
#' The return values of the agglomeration function are returned for each value 
#' of \eqn{i} in a data frame.
#' @return Returns a data frame with one line for every \eqn{i} in the set 
#' specified by \code{subtrack.length}. The first column contains the values 
#' of \eqn{i}, the following columns contain the values of the agglomeration 
#' function applied on the measure values of tracks of exactly \eqn{i} segments.
#' @examples 
#' require(ggplot2)
#' dat <- aggregate(TCells, displacement, FUN="mean.se")
#' ggplot(dat, aes(x=i, y=mean)) + geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)
aggregate.tracks <- function( x, measure, by="subtracks", FUN=mean, 
    subtrack.length=seq(1, (maxTrackLength(x)-1)),
    max.overlap=max(subtrack.length), na.rm=FALSE, ... ){
    if( class( x ) != "tracks" ){
    	if( class( x ) %in% c("data.frame","matrix" ) ){
    		tracks <- wrapTrack( x )
    	} else {
    		stop("Cannot coerce argument 'tracks' to tracks object" )
    	}
    }
    if (max(subtrack.length) > (maxTrackLength(x)-1)) {
		warning("No track is long enough!")
  	}
  	if (is.character(FUN)) {
		if (FUN == "mean.ci.95") {
		  the.statistic <- function(x,...) {
		  	mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x)$conf.int, error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mx, upper=ci[2]))
		  }
		} else if (FUN == "mean.ci.99") {
		  the.statistic <- function(x,...) {
			mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x, conf.level=.99)$conf.int, 
				error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
		  }
		} else if (FUN == "mean.se") {
		  the.statistic <- function(x,...) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)/sqrt(sum(!is.na(x))) 
		  		# note that this also works is na.omit is FALSE
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "mean.sd") {
		  the.statistic <- function(x,...) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "iqr") {
		  the.statistic <- function(x,...) {
			lx <- quantile(x,probs=c(.25,.5,.75),na.rm=na.rm)
			names(lx) <- c("lower","median","upper")
			return(lx)
		  }
		} else {
			the.statistic <- match.fun( FUN )
		}
	} else {
		if ( !is.function( FUN) ) {
		   stop("FUN must be a function or a string")
		} else {
			the.statistic <- FUN
		}
	}
	if( by == "subtracks" ){
		# optimized version: avoids making copies of subtracks (relevant especially
		# for measures like displacement, which actually only consider the track
		# endpoints)
		if( length( intersect(c("track","limits"),names(formals(measure))) ) == 2 ){
			measure.values <- lapply( subtrack.length, 
				function(i) .ulapply( Filter(function(t) nrow(t)>i, x),
					function(t) apply( .subtrackIndices(t,i,min(max.overlap,i-1)), 
						 1, measure, track=t ) ) )
		} else {
		# unoptimized version: makes copies of all subtracks.
			the.subtracks <- list()
			measure.values <- list()
			for (i in subtrack.length) {
				the.subtracks[[i]] <- subtracks(x, i, min(max.overlap,i-1) ) 
				the.subtracks[[i]] <- Filter( Negate(is.null), the.subtracks[[i]] )
				measure.values[[i]] <- sapply( the.subtracks[[i]], measure )
			}
		}
	} else if (by == "prefixes" ) {
		measure.values <- lapply( subtrack.length, 
				function(i) sapply( Filter( Negate(is.null), prefixes( x, i )  ), 
					measure ) 
			)
	} else {
		stop( "grouping method ", by, " not supported" )
	}
	measure.values <- Filter( Negate(is.null), measure.values)
	value <- sapply(measure.values, the.statistic)
	ret <- rbind(i=subtrack.length, value)
	return(data.frame(t(ret)))
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

# TODO include ellipse
# TODO refer to PNAS paper
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
#' @param step.spacing How many positions are to be left out between 
#' the steps that are considered for the test. For persistent motion, subsequent
#' steps will be correlated, which leads to too low p-values because Hotelling's 
#' test assumes that the input data is independent. To avoid this, the resulting
#' p-value should either be corrected for this dependence (e.g. by adjusting 
#' the degrees of freedom accordingly), or `step.spacing` should be set to a value
#' high enough to ensure that the considered steps are approximately independent.
#' @param ... further arguments passed on to \code{plot}.
#' 
#' @details Computes the displacement vectors of all segments in the tracks 
#' given in \code{tracks}, and performs Hotelling's T-square Test on that vector.
#' (see \code{\link[DescTools]{HotellingsT2Test}}).
#' 
#' @examples 
#' ## Test H_0: T-cells are migrating by random walk on x and y coordinates,
#' ## and report the p-value
#' 
#' if( require( DescTools ) ){
#'    hotellingsTest( TCells, plot=FALSE )$p.value
#' }

hotellingsTest <- function(tracks, dim=c("x", "y"), 
	step.spacing=0, plot=FALSE, ... ) {
  if( !requireNamespace("DescTools",quietly=TRUE) ){
    stop("This function requires the 'DescTools' package.")
  }
  if( !requireNamespace("DescTools",quietly=TRUE) ){
    stop("This function requires the package 'DescTools'.")
  }
  stopifnot( !plot || (length(dim) %in% c(1,2) ))
  Tx <- projectDimensions(tracks, dim)
  sx <- t(sapply( subtracks(Tx, 1, overlap=-step.spacing), displacementVector ) )
  if (length(dim) > 1) {
    if (plot) {
      plot(sx,...)
      points(0, 0, col=3)
      points(mean(sx[, 1:ncol(sx)]), col=2, pch=4, cex=2)
    }
  } else {
    if (plot) {
      mean.sx <- mean(sx)
      sd.sx <- sd(sx)
      plot(0, 0, col=3, xlab = dim[1],...)
      stripchart(sx, add=TRUE, at=0, pch=1)
      points(mean.sx, 0, col=2, pch=4, cex=2)
      points(0, 0, col=3)
      print(mean.sx -  sd.sx)
      points(mean.sx - sd.sx / sqrt(length(sx)), 0, pch=3, col=4)
      points(mean.sx + sd.sx / sqrt(length(sx)), 0, pch=3, col=4)
    }
  }
  return(DescTools::HotellingsT2Test(sx))
}

#' Extract Tracks Based on User-Defined Properties
#' 
#' Given a tracks object, extract a subset based on upper and lower bounds of a certain
#' measure. For instance, extract all tracks with a certain minimum length.
#'
#' @param tracks the tracks from which subsetting is to be done on 
#' @param measure the track measure for which the subsetting is to be based on 
#' @param ll specifies the lower-limit of the allowable measure
#' @param ul specifies the upper-limit of the allowable measure
getSubsetOfTracks <- function(tracks,measure,ll,ul){
  if(class(tracks)=="tracks"){
    mdat <- sapply(tracks,measure)
    selection <- tracks[(mdat<=ul & mdat>=ll) == T]
    return(as.tracks(selection))
  }else{
    obj <-deparse(substitute(tracks))
    warning(sprintf("Obj:%s is not a tracks object",obj))
  }
}

#' Plot tracks in 3D
#'
#' Takes an input tracks object and plots them in 3D using the 
#' \link[scatterplot3d]{scatterplot3d} function.
#'
#' @param tracks the tracks which will be plotted in 3d
#' @param ... further arguments to be passed on to 
#' \link[scatterplot3d]{scatterplot3d}
#'
plot3d <- function(tracks,...){
	if( !requireNamespace("scatterplot3d",quietly=TRUE) ){
		stop("This function requires the package 'scatterplot3d'.")
	}
	tracks_df <- as.data.frame.tracks(
		lapply( tracks, function(t) rbind(t,rep(NA,ncol(t)) ) ) )
	s3d <- scatterplot3d::scatterplot3d(tracks_df[,-c(1,2)],
		type="n",xlab="X Position",ylab="Y Position",zlab="Z Position",...)
	colvec <- rainbow(length(names(tracks)))[tracks_df[,1]]
	pts <- s3d$xyz.convert( tracks_df[,-c(1,2)] )
	segments( head(pts$x,-1), head(pts$y,-1), tail(pts$x,-1), 
		tail(pts$y,-1), col=colvec )
}

#' Compute Time Step of Tracks
#'
#' Applies a summary statistics on the time intervals between pairs of consecutive 
#' positions in a track dataset. 
#' 
#' @param tracks the input tracks
#' @param FUN the summary statistic to be applied.
#'
#' @details Most track quantification depends on the assumption that track positions are
#' recorded equidistantly in time. If this is not the case, then most of the statistic
#' in this (except for some very simple ones like \code{\link{duration}}) will not work.
#' In reality, at least small fluctuations of the time steps can be expected. This
#' function provides a means for quality control with respect to the tracking time.
#'
#' @return summary statistic of the time between two consecutive positions in a track 
#' dataset. 
#'
#' @examples
#' ## Show tracking time fluctuations for the T cell data
#' d <- timeStep( TCells )
#' plot( sapply( subtracks( TCells, 1 ) , duration ) - d, ylim=c(-d,d) ) 
timeStep <- function(tracks, FUN=mean ){
	if( class(FUN) == "character" ){
		FUN=match.fun( FUN )
	}
	FUN( sapply( subtracks( tracks, 1 ) , duration )  )
}

