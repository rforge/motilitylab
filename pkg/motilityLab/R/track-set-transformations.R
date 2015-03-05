#' Project Tracks onto given Dimensions
#' 
#' Projects tracks in a \emph{tracks} object onto the given dimensions.
#' 
#' @param tracks the \emph{tracks} onject whose tracks are to be projected.
#' @param dims a vector of characters giving the dimensions onto which the 
#' tracks are to be projected.
#' 
#' @details Returns only the given dimensions of the input tracks.
#' 
#' @return A tracks object is retruned which contains only those dimensions 
#' of the input \code{tracks} that are given in \code{dims}.
projectDimensions <- function(tracks, dims=c("x","y")) {
	as.tracks(lapply(tracks, function(t) {
		t[, c("t",dims)]
	}))
}
