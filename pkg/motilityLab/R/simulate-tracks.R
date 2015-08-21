#' Simulate \code{nsteps} Steps of a Random Walk in \code{dim} Dimensions
#' 
#' Generates a random track with \code{nsteps} steps in \code{dim} dimensions.
#' 
#' @param nsteps number of steps the track shall consist of.
#' @param dim the number of dimensions the track shall have.
#' @param mean stepwise mean drift per dimension; use 0 for an 
#' unbiased Brownian motion and other values for Brownian motion with drift
#' @param sd stepwise standard deviation per dimension
#' 
#' @details In in every step an for each dimension, a normally distributed 
#' value with mean \code{mean} and standard deviation \code{sd} is 
#' added to the previous cell position.
#' 
#' @return A data frame  containing in cell track with \code{nsteps} steps in 
#' \code{dim} dimensions is returned.
brownianTrack <- function(nsteps=100, dim=3, mean=0, sd=1, ...) {
	m <- cbind( 0:nsteps, replicate(dim, diffinv(rnorm(nsteps,mean=mean,sd=sd))))
	colnames(m) <- c("t", "x", "y", "z")
	m
}

#' Generate Tracks by Simulation
#'
#' Generic function that executes \code{expr}, which is expected to 
#' return a track, \code{n} times and stores the output in a \code{tracks}
#' object. Basically, this works like \code{\link{replicate}} but for tracks.
#'
#' @param n number of tracks to be generated. 
#' @param expr the expression, usually a call, that generates a single track.
#'
#' @return a \code{tracks} object containing \code{n} tracks.
#'
#' @example
#' ## Generate 10 tracks, 100 steps each, from a random walk with standard normally
#' ## distributed increments and plot them
#' plot( simulateTracks( 10, brownianTrack(100,3) ) )
simulateTracks <- function( n, expr ){
	as.tracks( 
		sapply( as.character(seq_len(n)), eval.parent(substitute(function(...) expr)), 
		simplify=FALSE, USE.NAMES=TRUE ) )
}
