#' Simulate \code{nsteps} Steps of a Random Walk in \code{dim} Dimensions
#' 
#' Generates a random track with \code{nsteps} steps in \code{dim} dimensions.
#' 
#' @param nsteps number of steps the track shall consist of.
#' @param dim the number of dimensions the track shall have.
#' 
#' @details In in every step an for each dimension, a normally distributed 
#' value is added to the previous cell position.
#' 
#' @return A data frame  containing in cell track with \code{nsteps} steps in 
#' \code{dim} dimensions is returned.
brownianTrack <- function(nsteps=100, dim=3) {
	m <- cbind( 0:nsteps, replicate(dim, diffinv(rnorm(nsteps))))
	colnames(m) <- c("t", "x", "y", "z")
	as.data.frame(m)
}

