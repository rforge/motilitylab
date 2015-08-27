#' Extract Subset of Dimensions
#' 
#' Projects tracks in a \emph{tracks} object onto the given dimensions.
#' 
#' @param tracks the \emph{tracks} onject whose tracks are to be projected.
#' @param dims a character vector giving the dimensions to extract from each track.
#' 
#' @return A tracks object is returned that contains only those dimensions 
#' of the input \code{tracks} that are given in \code{dims}.
projectDimensions <- function(tracks, dims=c("x","y")) {
	as.tracks(lapply(tracks, function(t) {
		t[, c("t",dims)]
	}))
}
