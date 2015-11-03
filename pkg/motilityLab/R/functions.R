#' Convert Tracks to Data Frame
#' 
#' Converts tracks from the list-of-matrices format, which is good
#' for efficient processing and therefore the default in this package, to a 
#' single dataframe which is convenient for plotting or saving the data.
#' 
#' @param x the \code{tracks} object to be coerced to a data frame.
#' 
#' @param row.names NULL or a character vector giving row names for the
#'  data frame.  Missing values are not allowed.
#' @param optional logical. Required for S3 consistency, but 
#' has no effect: column names are always assigned to the resulting
#'  data frame regardless of the setting of this option.
#' @param include.timepoint.column logical. If set to \code{TRUE}, then the resulting
#'  dataframe will contain a column that consecutively numbers the positions according
#'  to their time. Note that this information is anyway implicitly present in the time 
#'  information.
#' @param ... further arguments to be passed from or to other methods.
#'
#' @return a single data frame containing all individual tracks from the input with a 
#' prepended column named "id" containing each track's identifier in `x`.
#'
#' @examples
#' ## Display overall average position of the T cell data
#' colMeans( as.data.frame( TCells )[-c(1,2)] )
as.data.frame.tracks <- function(x, row.names = NULL, optional = FALSE, 
	include.timepoint.column=FALSE, ...) {
	ids <- rep(names(x), lapply(x,nrow))
	if( include.timepoint.column ){
		timepoint <- ave( ids, ids, FUN=seq_along )
		r <- data.frame(id=ids, timepoint=timepoint, do.call( 
			rbind.data.frame, x ) )
	} else {
		r <- data.frame(id=ids, do.call( 
			rbind.data.frame, x ) )
	}
	if( !is.null(row.names) && length(row.names)==nrow(r) ){
		rownames(r) <- row.names
	} else {
		rownames(r) <- NULL
	}
	return(r)
}

`[.tracks` <- function(x,y) as.tracks(as.list(x)[y])

#' Convert to Tracks
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
#' "Converts" a list into a \code{tracks} object by settings its 
#' S3 class to \code{tracks}.
#' 
#' @param x the list that is to be converted into a \emph{tracks} object.
#' @param ... further arguments to be passed from or to other methods.
#' 
#' @return the coerced object of S3 class \code{tracks} (i.e., a list 
#'  containing a matrix for each track)
as.tracks.list <- function(x,...) 
  structure(x, class="tracks")

#' Check for Class Tracks
#'
#' Verifies whether the input object inherits from S3 class \code{"tracks"}.
#'
#' @param x the input object.
is.tracks <- function(x)
  inherits(x, "tracks")

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

#' Sort Track Positions by Time
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
	as.tracks(lapply(x, function(d) d[order(d[,"t"],decreasing=decreasing),,drop=FALSE] ))
}


#' Join Track Datasets
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


#' Read Tracks from Text File
#' 
#' Reads cell tracks from a CSV or other text file. Data are expected to be organized as
#' follows.
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
#' the spatial coordinates. If this vector has two entries and the second entry is NA,
#' e.g. \code{c('x',NA)} or \code{c(5,NA)} then all columns from the indicated column 
#' to the last column are used. This is useful when reading files where the exact number
#' of spatial dimensions is not known beforehand.
#'
#' @param scale.t a value by which to multiply each time point. Useful for changing units,
#' or for specifying the time between positions if this is not contained in the file
#' itself.
#'
#' @param scale.pos a value, or a vector of values, by which to multiply each spatial 
#' position. Useful for changing units.
#' 
#' @param header a logical value indicating whether the file contains the
#' names of the variables as its first line. See \code{\link[utils]{read.table}}.
#'
#' @param sep a character specifying how the colums of the data are separated. 
#' The default value \code{""} means columns are separated by tabs or other spaces.  
#' See \code{\link[utils]{read.table}}.
#'
#' @param track.sep.blankline logical. If set to \code{TRUE}, then tracks are expected
#' to be separated by one or more blank lines in the input file instead of being 
#' designated by a track ID column. In this case, numerical track IDs are automatically
#' generated.
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
	pos.columns=c(3,4,5), scale.t=1, scale.pos=1, 
	header=TRUE, sep="", track.sep.blankline=FALSE, ...) {
	if( track.sep.blankline ){
		data.raw <- read.table(file, header=header, sep=sep,
			blank.lines.skip = FALSE, ...)
		is.blank <- apply( data.raw, 1, function(x) all( is.na(x) ) )
		ids <- cumsum( is.blank )
	} else {
		data.raw <- read.table(file, header=header, sep=sep, ...)
	}
	if( ncol(data.raw) < length(pos.columns)+1+as.integer(!track.sep.blankline) ){
		stop("CSV file does not contain enough columns! (Perhaps you need to specify 'sep')")
	}
	if( length(pos.columns) < 1 ){
		stop("At least one position column needs to be specified!")
	}
	if( length(pos.columns) == 2 && !is.finite(pos.columns[2]) ){
		cx <- match( pos.columns[1], colnames( data.raw ) )
		if( is.na(cx) && is.numeric(pos.columns[1]) ){
			cx <- pos.columns[1]
		}
		pos.columns <- seq( cx, length(data.raw) )
	}
	if( track.sep.blankline ){
		data.raw <- cbind( data.raw, ids )
		data.raw <- data.raw[!is.blank,]
		id.column <- ncol( data.raw )
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

#' Split Track into Multiple Tracks
#' 
#' @param x the input track (a data frame or a matrix).
#' @param id a string used to identify the resulting tracks; for instance,
#' if \code{id="X"}, then the resulting tracks are named X_1, X_2 and so forth.
#' Otherwise, they are simply labelled with integer numbers.
#' @param positions a vector of positive integers, given in ascending order.
#' @param min.length nonnegative integer. Resulting tracks that have fewer positions than
#'  the value of this parameter are dropped. 
splitTrack <- function( x, positions, id=NULL, min.length=2 ){
	freqs <- diff(c(0,positions,nrow(x)))
	segs <- rep( seq_len(length(positions)+1), freqs )
	segs[rep(freqs,freqs) < min.length] <- NA
	r <- structure( 
		split.data.frame( as.matrix(x), segs ),
		class="tracks")
	if(!is.null(id) && ( length(r) > 0 ) ){
		if( length(r) > 1 ){
			names( r ) <- paste0( id, "_", seq_along(r) )
		} else if( length(r) == 1 ){
			names( r ) <- id
		}
	}
	return(r)
}

#' Plot Tracks in 2D
#' 
#' Plots tracks contained in a "tracks" object into a twodimensional space
#' pallelel to the data's axes.
#' 
#' @param x the tracks to be plotted.
#' @param dims a vector giving the dimensions of the track data that shall be 
#' plotted, e.g. \code{c('x','y')} for the \eqn{x} and \eqn{y} dimension.
#' @param add boolean value indicating whether the tracks are to be added to the
#' current plot.
#' @param col a specification of the color(s) to be used. This can be a vector
#'  of size \code{length(x)}, where each entry specififes the color for the 
#'  corresponding track.
#' @param pch.start point symbol with which to label the first position of the track
#'  (see \code{\link[graphics]{points}}).
#' @param pch.end point symbol with which to label the last position of the track
#' @param cex point size for positions on the tracks.
#' @param ... additional parameters (e.g. xlab, ylab).
#' to be passed to \code{\link[graphics]{plot}}
#' (for \code{add=FALSE}) or \code{\link[graphics]{points}} (for \code{add=TRUE}), 
#' respectively.
#' @details One dimension of the data (by default \eqn{y}) is plotted against 
#' another (by default \eqn{x}). The dimesions can be chosen by means of the 
#' parameter \code{dims} and the axes can be labeled accordingly with the aid 
#' of \code{xlab} and \code{ylab}. The color can be set through \code{col}.
#' If the tracks should be added to an existing plot, \code{add} is to be set 
#' to \code{TRUE}.
#' @seealso \code{\link{plot3d}}
plot.tracks <- function(x, dims=c('x','y'), add=F,
		col=order(names(x)), pch.start=1, pch.end=NULL, 
		cex=.5, ... ) {
	args <- list(...)

	ids <- names(x)
	lcol <- rep_len(col, length(ids))
	pcol <- rep(lcol, lapply(x,nrow))
	dim1 <- unlist( lapply(x,'[', , dims[1]) )
	dim2 <- unlist( lapply(x,'[', , dims[2]) )

	if (add==T) {
		points( dim1, dim2, col=pcol, cex=cex, ... )
	} else {
		plot( dim1, dim2, col=pcol, cex=cex, ...)
	}
  
	if( !is.null(pch.start) ){
		starting.points <- lapply( x, function(p) p[1,] )
		px <- sapply(starting.points,'[[',dims[1])
		py <- sapply(starting.points,'[[',dims[2])
		points( px, py, col=col, pch=pch.start, cex=2 )
	}

	if( !is.null(pch.end) ){
		end.points <- lapply( x, function(p) p[nrow(p),] )
		px <- sapply(end.points,'[[',dims[1])
		py <- sapply(end.points,'[[',dims[2])
		points( px, py, col=col, pch=pch.end, cex=2 )
	}
  
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

#' Create Track Dataset from Single Track
#' 
#' Makes a \code{tracks} object containing the given track.
#' 
#' @param x the input track.
#' 
#' @return a list of class \code{tracks} containing only the input track \code{x}, which
#' is assigned the name "1". 
wrapTrack <- function(x) {
  return(as.tracks(list("1" = x)))
}

#' Staggered Version of a Function
#' 
#' Returns the "staggered" version of a track measure. That is, instead of
#' computing the measure on the whole track, the measure is averaged over
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
#' @references 
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
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
#' ## Compute the staggered matrix for overallAngle applied to all long enough 
#' ## subtracks of the first T cell track
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


#' Measure for Every Prefix of a Track
#' 
#' @param measure the measure that is to be computed.
#' @param min.length minimal number of points in a prefix the measure is to be 
#' applied to. Default is \eqn{1}.
#' 
#' @details The resulting 
#' function applies \code{measure} to every prefix of the input track, starting
#' with the prefix that consists of \code{min.length} positions. The resulting 
#' values are returned in a vector, sorted in ascending order by the 
#' prefix length.
#'
#' @return a function that applies the measure to every prefix of length at 
#' least \code{min.length}.
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

#' Decompose Track(s) into Subtracks
#' 
#' Creates a \emph{tracks} object consisting of all subtracks of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \code{tracks} object.
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
	if( !is.tracks(x) ){
		x <- wrapTrack( x )
	}
	do.call(c, lapply(x, 
  		function(t) .computeSubtracksOfISegments(t, i, overlap )))
}

#' Get Track Prefixes 
#' 
#' Creates a \code{tracks} object consisting of all prefixes (i.e., subtracks
#' starting with the first position of a track) of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \code{tracks} object.
#' @param i subtrack length. A single integer, lists are not supported.
#'
#' @details This function behaves exactly like \code{\link{subtracks}} except
#' that only subtracks starting from the first position are considered.
#'
prefixes <- function(x,i) {
	if( !is.tracks(x) ){
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
#' @return a function that computes the input measure for a given track
#' and returns the result divided by the track's duration.
#' 
#' @examples 
#' ## normalizeToDuration(displacement) can be used as an indicator 
#' ## for the motion's efficiency
#' sapply(TCells, normalizeToDuration(displacement))
normalizeToDuration <- function(measure) {
  function(track) {
    measure(track) / duration(track)
  }
}



#' Bivariate Scatterplot of Two Given Track Measures
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
#' applied to using \code{sapply} before plotting the result.
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
#' value of  \code{measurey} applied to \code{y}. This is useful for "FACS-like"
#' motility analysis, where clusters of cell tracks are identified based on their
#' motility parameters (Moreau et al, 2012; Textor et al, 2014).
#'
#' @references
#' Moreau HD, Lemaitre F, Terriac E, Azar G, Piel M, Lennon-Dumenil AM,
#' Bousso P (2012), Dynamic In Situ Cytometry Uncovers 
#' T Cell Receptor Signaling during Immunological Synapses and Kinapses In Vivo.
#' \emph{Immunity} \bold{37}(2), 351--363. doi:10.1016/j.immuni.2012.05.014
#'
#' Johannes Textor, Sarah E. Henrickson, Judith N. Mandl, Ulrich H. von Andrian,
#' J\"urgen Westermann, Rob J. de Boer and Joost B. Beltman (2014),
#' Random Migration and Signal Integration Promote Rapid and Robust T Cell Recruitment.
#' \emph{PLoS Computational Biology} \bold{10}(8), e1003752.
#' doi:10.1371/journal.pcbi.1003752
#'
plotTrackMeasures <- function(x, measurex, measurey, add=FALSE, 
                              xlab=deparse(substitute(measurex)), 
                              ylab=deparse(substitute(measurey)), ...) {
	if( !is.tracks(x) ){
		x <- wrapTrack( x )
	}
    if(add) {
      points(sapply(y,measurey) ~ sapply(x,measurey), ...)
    } else {
      plot(sapply(y,measurey) ~ sapply(x,measurey), xlab=xlab, ylab=ylab, ...)
    }
}
  

#' Cluster Tracks
#' 
#' Perform a hierarchical clustering of a set of tracks according to a given vector
#' of track measures.
#' 
#' @param tracks the tracks that are to be clustered.
#' @param measures a function, or a vector of functions. Each function is expected to 
#' return a single number given a single track.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the clustering. 
#' @param ... additional parameters to be passed to \code{\link[stats]{hclust}}.
#'
#' @return An object of class *hclust*, see \code{\link[stats]{hclust}}.
#' 
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object. According to the resulting values, the tracks are 
#' clustered using a hierarchical clustering (see \code{\link[stats]{hclust}}).
#' If \code{scale} is \code{TRUE}, the measure values are scaled to mean value 
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the clustering.
#'
#' @examples 
#' ## Cluster tracks according to the mean of their Hust exponents along X and Y
#'
#' cells <- c(TCells,Neutrophils)
#' real.celltype <- rep(c("T","N"),c(length(TCells),length(Neutrophils)))
#' ## Prefix each track ID with its cell class to evaluate the clustering visually
#' names(cells) <- paste0(real.celltype,seq_along(cells))
#' clust <- clusterTracks( cells, hurstExponent )
#' plot( clust )
#' ## How many cells are "correctly" clustered?
#' sum( real.celltype == c("T","N")[cutree(clust,2)] )
#' 
clusterTracks <- function(tracks, measures, scale=TRUE, ...) { 
	values <- matrix(nrow=length(tracks))
	if (is.function(measures)) {
		measures <- c(measures)
	}
	values <- do.call(cbind, lapply(measures, function(m) sapply(tracks, m)))
	if (scale) {
		values <- scale(values)
	}
	hclust(dist(values), ...)
}


#' Compute Summary Statistics of Subtracks
#' 
#' Computes a given measure on subtracks of a given track set, applies a summary 
#' statistic for each subtrack length, and returns the results in a convenient form.
#' This is the main workhorse function that facilitates the most common motility analyses 
#' such as mean square displacement, turning angle, and autocorrelation plots.
#' 
#' @param x the tracks object whose subtracks are to be considered. 
#' If a single track is given, it will be coerced to a tracks object 
#' using \code{\link{wrapTrack}} (but note that this requires an explicit call
#' \code{aggregate.tracks}).
#' @param measure the measure that is to be computed on the subtracks.
#' @param by a string that indicates how grouping is performed. Currently, two
#' kinds of grouping are supported: 
#' \itemize{
#'  \item{"subtracks"}{ Apply \code{measure} to all subtracks according to 
#' the parameters \code{subtrack.length} and \code{max.overlap}.}
#'  \item{"prefixes"}{ Apply \code{measure} to all prefixes (i.e., subtracks starting
#'  from a track's initial position) according to the parameter \code{subtrack.length}.}
#' }
#' 
#' @param FUN a summary statistic to be computed on the measures of subtracks with the
#' same length. Can be a function or a string.
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
#' \code{\link{match.fun}}.
#'
#' @param subtrack.length an integer or a vector of integers defining which subtrack
#' lengths are considered. In particular, \code{subtrack.length=2} 
#' corresponds to a "step-based analysis" (Beltman et al, 2009).
#'
#' @param max.overlap an integer controlling what to do with overlapping subtracks.
#' A maximum overlap of \code{max(subtrack.length)} will imply
#' that all subtracks are considered. For a maximum overlap of 0, only non-overlapping
#' subtracks are considered. A negative overlap can be used to ensure that only subtracks
#' a certain distance apart are considered. In general, for non-Brownian motion there will
#' be correlations between subsequent steps, such that a negative overlap may be necessary
#' to get a proper error estimate.
#' @param na.rm logical. If \code{TRUE}, then \code{NA}'s and \code{NaN}'s 
#' are removed prior to computing the summary statistic.
#' @param filter.subtracks a function that can be supplied to exclude certain subtracks
#' from an analysis. For instance, one may wish to compute angles only between steps of
#' a certain minimum length (see Examples).
#' @param ... further arguments passed to or used by methods.
#'
#' @details For every number of segments \eqn{i} in the set defined by 
#' \code{subtrack.length}, all subtracks of any track in the input 
#' \code{tracks} object that consist of exactly \eqn{i} segments are 
#' considered. The input \code{measure} is applied to the subtracks individually, 
#' and the \code{statistic} is applied to the resulting values. 
#' 
#' @return a data frame with one row for every \eqn{i} 
#' specified by \code{subtrack.length}. The first column contains the values 
#' of \eqn{i} and the remaining columns contain the values of the summary statistic
#' of the measure values of tracks having exactly \eqn{i} segments.
#' @examples 
#' ## A mean square displacement plot with error bars.
#' dat <- aggregate(TCells, squareDisplacement, FUN="mean.se")
#' with( dat ,{
#'   plot( mean ~ i, xlab="time step", 
#'   	ylab="mean square displacement", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
#'
#' ## Compute a turning angle plot for the B cell data, taking only steps of at least
#' ## 1 micrometer length into account
#' check <- function(x) all( sapply( list(head(x,2),tail(x,2)), trackLength ) >= 1.0 )
#' plot( aggregate( BCells, overallAngle, subtrack.length=1:10, 
#'   filter.subtracks=check )[,2], type='l' )
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
aggregate.tracks <- function( x, measure, by="subtracks", FUN=mean, 
    subtrack.length=seq(1, (maxTrackLength(x)-1)),
    max.overlap=max(subtrack.length), na.rm=FALSE, 
    filter.subtracks=NULL, ... ){
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
  	the.statistic <- NULL
  	if (is.character(FUN)) {
		if (FUN == "mean.ci.95") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x)$conf.int, error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mx, upper=ci[2]))
		  }
		} else if (FUN == "mean.ci.99") {
		  the.statistic <- function(x) {
			mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x, conf.level=.99)$conf.int, 
				error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
		  }
		} else if (FUN == "mean.se") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)/sqrt(sum(!is.na(x))) 
		  		# note that this also works is na.omit is FALSE
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "mean.sd") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "iqr") {
		  the.statistic <- function(x) {
			lx <- quantile(x,probs=c(.25,.5,.75),na.rm=na.rm)
			names(lx) <- c("lower","median","upper")
			return(lx)
		  }
		} else {
			the.statistic.raw <- match.fun( FUN )
		}
	} else {
		if ( !is.function( FUN) ) {
		   stop("FUN must be a function or a string")
		} else {
			the.statistic.raw <- FUN
		}
	}
	if( is.null( the.statistic ) ){
		if( na.rm ){
			the.statistic <- function(x){
				the.statistic.raw(x[!is.na(x) & !is.nan(x)])
			}
		} else {
			the.statistic <- the.statistic.raw
		}
	}
	if( is.function( filter.subtracks ) ){
		isValidSubtrack <- function(t) !is.null(t) && filter.subtracks(t)
	} else {
		isValidSubtrack <- Negate(is.null)
	}
	if( by == "subtracks" ){
		# optimized version: avoids making copies of subtracks (relevant especially
		# for measures like displacement, which actually only consider the track
		# endpoints). Also pre-computes the diff of the track if needed.
		if( ( length( intersect(c("x","limits"),names(formals(measure))) ) == 2 )
			&& !is.function( filter.subtracks ) ){
			if( "xdiff" %in% names(formals(measure)) ){
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)), 
					1, measure, x=t, xdiff=diff(t) )
			} else {
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)), 
					1, measure, x=t )
			}
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				measure.values[[k]] <- .ulapply( xi, ft, i=i )
				k <- k+1
			}
		} else if( ( length( intersect(c("x","from","to"),names(formals(measure))) ) == 3 )
			&& !is.function( filter.subtracks ) ){
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				xi1 <- do.call( rbind, xi )
				si <- .subtrackIndices(i,min(max.overlap,i-1),sapply(xi,nrow))
				measure.values[[k]] <- measure(xi1,si[,1],si[,2])
				k <- k+1
			}
		} else {
			# unoptimized version: makes copies of all subtracks.
			the.subtracks <- list()
			measure.values <- list()
			k <- 1
			for (i in subtrack.length) {
				the.subtracks[[k]] <- subtracks(x, i, min(max.overlap,i-1) )
				the.subtracks[[k]] <- Filter( isValidSubtrack, the.subtracks[[k]] )
				measure.values[[k]] <- sapply( the.subtracks[[k]], measure )
				k <- k+1
			}
		}
	} else if (by == "prefixes" ) {
		measure.values <- lapply( subtrack.length, 
				function(i) sapply( Filter( isValidSubtrack, prefixes( x, i )  ), 
					measure ) 
			)
	} else {
		stop( "grouping method ", by, " not supported" )
	}
	# this is necessary because of automatic insertion of NULL elements when e.g.
	# measure.values[[5]] is assigned to
	measure.values <- Filter( Negate(is.null), measure.values )
	value <- sapply(measure.values, the.statistic)
	ret <- rbind(i=subtrack.length, value)
	return(data.frame(t(ret)))
}

#' Length of Longest Track
#' 
#' Determines the maximum number of positions over the tracks in \code{x}.
#' 
#' @param x the \code{tracks} object the tracks in which are to be considered.
#' 
#' @return the maximum number of rows of a track in \code{x}
maxTrackLength <- function(x) {
  return(max(sapply(x, nrow)))
}


#' Bounding Box of a Tracks Object
#' 
#' Computes the minimum and maximum coordinates per dimension for all positions in a
#' given list of tracks. Assumes that all tracks have the same number of dimensions. 
#' 
#' @param x the \code{tracks} object whose bounding box is to be computed.
#'  
#' @return Returns a matrix with two rows and \eqn{d} columns, where \eqn{d} is 
#' the number of dimensions of the tracks. The first row contains the minimum 
#' and the second row the maximum value of any track in the dimension given by 
#' the column.
boundingBox <- function(x) {
	if( !is.tracks(x) ){
		x <- wrapTrack(x) 
	}
	bounding.boxes <- lapply(x, .boundingBoxTrack) 
	empty <- mat.or.vec(nrow(bounding.boxes[[1]]), ncol(bounding.boxes[[1]]))
	colnames(empty) <- colnames(bounding.boxes[[1]])
	rownames(empty) <- rownames(bounding.boxes[[1]])
	empty[1,] <- Inf
	empty[2,] <- -Inf
	return(Reduce(.minMaxOfTwoMatrices, bounding.boxes, empty)) 
}

#' Test Unbiasedness of Motion
#' 
#' Test the null hypothesis that a given set of tracks originates from an uncorrelated
#' and unbiased type of motion (e.g., a random walk without drift). This is done by
#' testing whether the mean step vector is equal to the null vector. 
#' 
#' @param tracks the tracks whose biasedness is to be determined.
#' @param dim vector with the names of the track's 
#' dimensions that are to be considered. By default c("x", "y").
#' @param plot logical indicating whether the scatter of the step's directions, 
#' origin of ordinates (green circle) and the mean of the data points (green 
#' cross) are to be plotted. (In one dimension also the bounds of the 
#' condfidence interval are given.) Plot works only in one or two dimensions.
#' @param add whether to add the plot to the current plot (\code{TRUE}) or create a 
#" new one (\code{FALSE}).
#' @param step.spacing How many positions are to be left out between 
#' the steps that are considered for the test. For persistent motion, subsequent
#' steps will be correlated, which leads to too low p-values because Hotelling's 
#' test assumes that the input data is independent. To avoid this, the resulting
#' p-value should either be corrected for this dependence (e.g. by adjusting 
#' the degrees of freedom accordingly), or `step.spacing` should be set to a value
#' high enough to ensure that the considered steps are approximately independent.
#' @param ellipse.col color with which to draw the confidence ellipse. Use \code{NA} to
#'  omit the confidence ellipse.
#' @param conf.level the desired confidence level for the confidence ellipse.
#' @param ... further arguments passed on to \code{plot}.
#' 
#' @return a list with class \code{htest}, see \code{\link[DescTools]{HotellingsT2Test}}.
#' 
#' @details Computes the displacement vectors of all segments in the tracks 
#' given in \code{tracks}, and performs Hotelling's T-square Test on that vector.
#' (see \code{\link[DescTools]{HotellingsT2Test}}).
#' 
#' @examples 
#' ## Test H_0: T-cells migrate by uncorrelated random walk on x and y coordinates,
#' ## and report the p-value.
#' if( require( DescTools ) ){
#'    hotellingsTest( TCells, plot=FALSE )$p.value
#' }
#'
#' @references
#' Johannes Textor, Antonio Peixoto, Sarah E. Henrickson, Mathieu
#'  Sinn, Ulrich H. von Andrian and Juergen Westermann (2011),
#'	Defining the Quantitative Limits of Intravital Two-Photon Lymphocyte Tracking.
#' \emph{PNAS} \bold{108}(30):12401--12406. doi:10.1073/pnas.1102288108
#'
hotellingsTest <- function(tracks, dim=c("x", "y"), 
	step.spacing=0, plot=FALSE, add=FALSE, ellipse.col="blue", conf.level=0.95, ... ) {
	if( !requireNamespace("DescTools",quietly=TRUE) ){
		stop("This function requires the package 'DescTools'. Please install it.")
	}
	stopifnot( !plot || (length(dim) %in% c(1,2) ))
	Tx <- projectDimensions(tracks, dim)
	sx <- t(sapply( subtracks(Tx, 1, overlap=-step.spacing), displacementVector ) )
	if (plot && (length(dim)==2)) {
		if (plot) {
			if( add ){
				points(sx,...)
			} else {
				plot(sx,...)
			}
			abline(h=0)
			abline(v=0)
			cm <- colMeans(sx)
			if( !is.na(ellipse.col) ){
				if( !requireNamespace("ellipse",quietly=TRUE) ){
					warning("Drawing confidence ellipse requires the package 'ellipse'. Please install it.")
				} else {
					n <- dim(sx)[1]
					p <- dim(sx)[2]
					S <- cov(sx)
					t <- sqrt(((n-1)*p/(n*(n-p)))*qf(1-conf.level,p,n-p,lower.tail=F))
					polygon(ellipse::ellipse(S,centre=cm,t=t),
						col=.setColAlpha(ellipse.col,128),border=NA)
					points(cm[1],cm[2],col=ellipse.col,pch=20)
				}
			}
		}
	} else if (plot && (length(dim)==1)) {
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
	return(DescTools::HotellingsT2Test(sx))
}

#' Select Tracks by Measure Values
#' 
#' Given a tracks object, extract a subset based on upper and lower bounds of a certain
#' measure. For instance, extract all tracks with a certain minimum length.
#'
#' @param x the input tracks.
#' @param measure measure on which the selection is based.
#' @param lower specifies the lower bound (inclusive) of the allowable measure.
#' @param upper specifies the upper bound (inclusive) of the allowable measure.
#' @examples
#' ## Slower half of T cells
#' slow.tcells <- selectTracks( TCells, speed, -Inf, median( sapply(TCells,speed) ) )
selectTracks <- function(x,measure,lower,upper){
	if(is.tracks(x)){
		mdat <- sapply(x,measure)
		return(x[mdat<=upper & mdat>=lower])
	}else{
		obj <-deparse(substitute(tracks))
		warning(sprintf("Obj:%s is not a tracks object",obj))
	}
}

#' Plot Tracks in 3D
#'
#' Takes an input tracks object and plots them in 3D using the 
#' \link[scatterplot3d]{scatterplot3d} function.
#'
#' @param x the tracks which will be plotted in 3d
#' @param ... further arguments to be passed on to 
#' \link[scatterplot3d]{scatterplot3d}
#'
#' @examples
#' if( require("scatterplot3d",quietly=TRUE) ){
#'   plot3d( TCells )
#' }
plot3d <- function(x,...){
	if( !requireNamespace("scatterplot3d",quietly=TRUE) ){
		stop("This function requires the package 'scatterplot3d'.")
	}
	tracks_df <- as.data.frame.tracks(
		lapply( x, function(t) rbind(t,rep(NA,ncol(t)) ) ) )
	s3d <- scatterplot3d::scatterplot3d(tracks_df[,-c(1,2)],
		type="n",xlab="X Position",ylab="Y Position",zlab="Z Position",...)
	colvec <- rainbow(length(names(x)))[tracks_df[,1]]
	pts <- s3d$xyz.convert( tracks_df[,-c(1,2)] )
	segments( head(pts$x,-1), head(pts$y,-1), tail(pts$x,-1), 
		tail(pts$y,-1), col=colvec )
}

#' Compute Time Step of Tracks
#'
#' Applies a summary statistics on the time intervals between pairs of consecutive 
#' positions in a track dataset. 
#' 
#' @param x the input tracks.
#' @param FUN the summary statistic to be applied.
#' @param na.rm logical, indicates 
#'  whether to remove missing values before
#'  applying FUN.
#'
#' @details Most track quantification depends on the assumption that track positions are
#' recorded at constant time intervals. If this is not the case, then most of the statistics
#' in this package (except for some very simple ones like \code{\link{duration}}) will not work.
#' In reality, at least small fluctuations of the time steps can be expected. This
#' function provides a means for quality control with respect to the tracking time.
#'
#' @return summary statistic of the time intervals between two consecutive positions in a track 
#' dataset. 
#'
#' @examples
#' ## Show tracking time fluctuations for the T cell data
#' d <- timeStep( TCells )
#' plot( sapply( subtracks( TCells, 1 ) , duration ) - d, ylim=c(-d,d) ) 
timeStep <- function( x, FUN=median, na.rm=FALSE ){
	if( class(FUN) == "character" ){
		FUN=match.fun( FUN )
	}
	x <- .ulapply( x, function(x) diff(x[,1]) )
	if( na.rm ){
		x <- x[!is.na(x)]
	}
	return( FUN( x ) )
}

#' Interpolate Track Positions
#'
#' Approximates the track positions at given time points using linear interpolation
#' (via the \code{\link[stats]{approx}} function).
#'
#' @param x the input track (a matrix or data frame).
#' @param t the times at which to approximate track positions. These must lie 
#' within the interval spanned by the track timepoints.
#' @param how specifies how to perform the interpolation. Possible values are 
#' \code{"linear"} (which uses \code{\link[stats]{approx}} with default values) and
#' \code{"spline"} (which uses \code{\link[stats]{spline}} with default values). 
#'
interpolateTrack <- function( x, t, how="linear" ){
	f=approx
	if( how=="spline" ){
		f=spline
	}
	pos.interpolated <- apply(x[,-1],2,function(p) f(x[,1], p, xout=t)$y) 
	r <- cbind( t, pos.interpolated )
	colnames(r) <- colnames(x)
	return(r)
}

#' Subsample Track by Constant Factor
#'
#' Make tracks more coarse-grained by keeping only every \emph{k}th position.
#'
#' @param x an input track or tracks object.
#' @param k a positive integer. Every \eqn{k}th position of each
#' input track is kept.
#' @seealso \code{interpolateTrack}, which can be used for more flexible track 
#' coarse-graining. 
#' @examples
#' ## Compare original and subsampled versions of the T cell tracks
#' plot( TCells, col=1 )
#' plot( subsample( TCells, 3 ), col=2, add=TRUE, pch.start=NULL )
subsample <- function( x, k=2 ){
	if( is.tracks(x) ){
		as.tracks( lapply( x, function(t) subsample(t,k) ) )
	} else {
		x[seq(1,nrow(x),by=k),]
	}
}
