#' Helper function of beauWalker() and is not directly called by user
#' This function samples uniformly from unit sphere
.beaucheminSphereSampler <- function(d=3){
  x <- rnorm(d)
  return(x/sqrt(sum(x^2)))
}

#' Helper function of beauWalker() and is not directly called by user
#' This function samples from the trianglular distribution
.beaucheminTriangularSampler <- function() sqrt(4*runif(1))-1

# Creates a matrix that will rotate unit vector a onto unit vector b
.beaucheminRotationMatrix <- function(a,b){
		# handle undefined cases
	    theta <- acos( sum(a*b) )
    	R <- diag(rep(1,3))
		if( theta < 0.001 ){
			return(R)
		}
		if( pi-theta < 0.001 ){
			return(-R)
		}
		# compute normalized cross product of a and b
		x <- c(a[2]*b[3]-a[3]*b[2],
			a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
		x <- x / sqrt(sum(x^2))
		A <- rbind( c(0,-x[3],x[2]),c(x[3],0,-x[1]),c(-x[2],x[1],0) )
		R <- R + sin(theta)*A + (1-cos(theta))* (A%*%A)
    	return(R)
}

#' Helper function of beauWalker() and is not directly called by user
#' returns a direction for the cell to travel based on model parameters specified in beauWalker()
.beaucheminPickNewDirection <- function(
	old.direction,
	p.bias,p.persist,bias.dir,taxis.mode,
	t.free,v.free,rot.mat){
	if( runif(1) < p.persist ){
		return(c(old.direction, t.free))
	}
	d <- .beaucheminSphereSampler(3)
	if( taxis.mode == 0 ){
		return(c(d,t.free))
	} else if(taxis.mode==1){
		# orthotaxis
		return(c(v.free * d*(1+p.bias*sum(d*bias.dir)),t.free))
	} else if(taxis.mode == 2){
		# topotaxis	
		if( runif(1) < p.bias ){
			# Approach: generate new direction as if the bias direction were (1,0,0).
			# Then rotate the resulting direction by the angles between (1,0,0) and the true
			# bias direction. 
			d[1] <- .beaucheminTriangularSampler()
			circ <- .beaucheminSphereSampler(2)
			circle.scale <- sqrt(1-d[1]^2)
			d[2:3] <- circ[1:2]*circle.scale
			return(c( v.free * rot.mat%*%d, t.free))
		}
	} else if(taxis.mode == 3){
	 	 # klinotaxis
		return(c(d*v.free,t.free*(1+p.bias*sum(d*bias.dir))))
	}
}

#' Simulates a single cell track based on updated Beauchemin model written by J. Textor, M. Sinn, and R. Boer. 
#' 
#' @param sim.time Specifies the duration time of the simulation (in minutes). When at time point i, if (i+1)*delta.t > sim.time
#' then the simulation ends. In otherwords, the actual simulation can only be either as long or just short of the allowed
#' simulation time. 
#' @param delta.t Change in time between each timepoint (in minutes).
#' @param p.persist Indicates how probable a change in direction is. p.persist = 1 indicates there is never a 
#' change in direction between steps and p.persist = 0 indicates there is always a change in direction between steps.
#' @param p.bias Strength of movement in the direction of \code{bias.dir}.
#' @param taxis.mode Specified mode of movement. 1 := orthotaxis, 2 := topotaxis, 3 := klinotaxis.
#' @param t.free Time interval for how long the cell is allowed to move between steps (in minutes).
#' @param v.free Speed of cell.
#' @param t.pause Time for how long the cell takes to adjust movement to new direction (in minutes).
#' @param bias.dir The three dimensional vector c(x,y,z) indicating the direction for which there is a preference for movement.
#' @details The default parameters are those found by mean-square displacement plots created by Beauchemin et al.
#' @return A matrix of size \code{n.steps}x5. Each row represents the three dimensional coordinates of the cell
#' at the given timepoint.
#' @examples 
#' ## Create track with model parameters and return matrix of positions
#' out <- beauWalker(n.steps=20,p.persist = 0.3,taxis.mode = 1)
#' ## Plot X-Y projection
#' plot(y~x,data=out,type='b',main="Planar projection over z-axis")
#' 
#' ## Create list of 20 tracks and plot all 
#' blist <- list()
#' for(i in 1:20){
#'  blist[[i]]<-beauWalker(bias.dir=c(-1,1,0),p.bias=10,taxis.mode = 1,p.persist = 0.1,n.steps = 10,delta.t = 1)
#" }
#' plot(x=NULL,xlim=c(-700,100),ylim=c(-100,700),type='b')
#" lapply(blist, function(i) points(y~x,data=i,type='b'))


beaucheminTrack <- function(sim.time=10,delta.t=1,p.persist=0.5,p.bias=0.9,bias.dir=c(0,0,0),taxis.mode=1,t.free=2,v.free=18.8,t.pause=0.5){  
  # Parameter checks
  if(p.persist < 0 || p.persist > 1){
    stop("p.persist must be a value between 0 and 1")
  }
  if(!(taxis.mode %in% 0:4)){
    stop("taxis.mode can either be 1, 2, 3, or 0 for no taxis. ",
    	"See documentation for details on which model you would like to use.")
  }
  if(p.bias < 0){
    stop("p.bias must be a value greater than or equal to 0")
  }
  if(delta.t <= 0){
    stop("delta.t must be a value larger than 0")
  }
  if(t.pause <= 0){
    stop("t.pause must be a value larger than 0")
  }
  if(t.free <= 0){
    stop("t.free must be a value larger than 0")
  }
  if(v.free <= 0){
    stop("v.free must be a value larger than 0")
  }

  # Initializing parameters 
  if(any(bias.dir != 0)){
	bias.dir <- bias.dir/sqrt(sum(bias.dir^2))
  }
  
  rot.mat <- NULL
  if(taxis.mode==2){
  	# cache rotation matrix for topotaxis to avoid recomputing it frequently
  	rot.mat <- .beaucheminRotationMatrix( c(1,0,0), bias.dir )
  }
  pos <- matrix(rep(0,4),1,4)
  
  n.steps <- ceiling( sim.time / (t.free+t.pause) ) + 1

  d <- .beaucheminPickNewDirection( NULL,p.bias,0,bias.dir,taxis.mode,
			t.free,v.free,rot.mat )
  p <- c(0,0,0)
  t <- t.pause
  for( i in seq_len(n.steps) ){
	pnew <- p+d[4]*d[-4]
	tnew <- t+d[4]
	pos <- rbind( pos, c(t,p), c(tnew,pnew) )
	t <- tnew + t.pause
	p <- pnew
	d <- .beaucheminPickNewDirection( d,p.bias,0,bias.dir,taxis.mode,
				t.free,v.free,rot.mat )
  }
  
  # interpolate track observations according to delta.t
  t <- seq(0,sim.time,by=delta.t)+runif(1,max=t.free+t.pause)
  pos.interpolated <- apply(pos[,-1],2,function(x) approx(pos[,1], x, xout=t)$y) 
  pos <- cbind( t, pos.interpolated )
  colnames(pos) <- c("t","x","y","z")

  return(normalizeTrack(pos))
}

#' A tracks fractal dimension
#' Computes the fractal dimension of a track using all track dimensions by the box-count method.
#' @param track the track for which the fractal dimension is to be calculated.
#' @details The fractal dimension is found using the function from \code{\link[fractaldim]{fd.estim.boxcount}} from the 
#' \code{fractaldim} package. While the hurst exponent takes a global approach to the track's properties, fractal dimension
#' is a local approach to the track's properties. These two are generally considered to be independent of one another.
#' @return the number indicating the track's fractal dimension.
fractalDimension <- function(track){
  if( !requireNamespace("fractaldim",quietly=TRUE) ){
    stop("This function requires the 'fractaldim' package.")
  }
  fd.full <- fractaldim::fd.estim.boxcount(track[,2:ncol(track)])
  fd.out <- fd.full$fd
  return(fd.out)
}




#' Removes constant drift from all tracks within a \code{tracks} object
#' 
#' @param x The \code{tracks} object for which the drift is to be removed from
#' @param drift A vector designating the direction of the drift
#' @details This function is to be used only if there is evidence that for each time point the movement is characterized 
#' by the equation:
#'    track[i+1,] = track[i,]+trueMovement+drift
#' Here drift is a constant vector that applies for all i thus the drift correction is simply subtracting out the designated 
#' drift vector from the positions of all time points excluding the initial time point.
#' @return A \code{tracks} object that has been corrected for drift
#' @examples 
#' blist <- list()
#' for(i in 1:20){
#'  blist[[i]]<-beauWalker(bias.dir=c(-1,1,0),p.bias=10,taxis.mode = 1,p.persist = 0.1,n.steps = 10,delta.t = 1)
#' }
#' blist.corrected <- correctDrift(blist,c(-1.5,.8,0))
correctDrift<- function(x, drift=NULL){
  if(is.null(drift)){
    stop("Please specify a vector indicating direction of the drift")
  }
  if(class(x) != "tracks"){
    stop("The supplied tracks argument is not a tracks object")
  }
  if(length(x)==0){
    stop("There are no tracks in the supplied input tracks object")
  }
  if(ncol(x[[1]])-1 != length(drift)){
    stop("The dimensions of drift do not match the dimensions of tracks")
  }
  return(as.tracks(lapply(x, function(i) rbind(i[1,],cbind(i[-1,1],i[2:nrow(i),2:ncol(i)]-drift)))))
}