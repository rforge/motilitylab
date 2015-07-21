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
	m
}


sbm.tracks <- function(ncells, randlength=FALSE, ...){

	tracks <- NULL

	if(randlength==FALSE){
 
		for (i in 1:ncells){
			tracks[[i]] <- brownianTrack(...)

		}

		names(tracks) <- as.character(seq(1:ncells))

		return(as.tracks(tracks))

	} else {

		for (i in 1:ncells){

			rnsteps <- floor(runif(n=1,min=1, max=101))

			tracks[[i]] <- brownianTrack(nsteps=rnsteps)

		}
 
		names(tracks) <- as.character(seq(1:ncells))
 
		return(as.tracks(tracks))

	}

}
	



gbm <- function(nsteps=100,dim=3,lag=0,drift=0){
  if(lag==0){
    for(i in 1:nsteps){
      xdis <- rnorm(nsteps, 0 ,1)
      ydis <- rnorm(nsteps, 0 ,1)
      zdis <- rnorm(nsteps, 0 ,1)
      xdis <- cumsum(xdis)+drift*nsteps
      ydis <- cumsum(ydis)+drift*nsteps
      zdis <- cumsum(zdis)+drift*nsteps
      m <- cbind(0:(nsteps-1),xdis,ydis,zdis)
    }
    colnames(m) <- c("t", "x", "y", "z")
    as.data.frame(m)
  }else{
    for(i in 1:nsteps){
      xdis <- rnorm(nsteps, 0 ,1)
      ydis <- rnorm(nsteps, 0 ,1)
      zdis <- rnorm(nsteps, 0 ,1)
      xdis <- cumsum(xdis)
      ydis <- cumsum(ydis)
      zdis <- cumsum(zdis)
      at <- rpois(nsteps,lag)
    }
    for(i in 1:nsteps){
      if(at[i]!=0){
      xdis[i] = xdis[i]*at[i]
      ydis[i] = ydis[i]*at[i]
      zdis[i] = zdis[i]*at[i]
      }
    }
    m <- cbind(0:(nsteps-1),xdis,ydis,zdis)
    colnames(m) <- c("t", "x", "y", "z")
    as.data.frame(m)
  }
}



gbm.tracks <- function(nsteps=100,ncells,randlength=F,lag=0,drift=0){
  tracks <- list()
  if(randlength==F){
    for (i in 1:ncells){
      tracks[[i]] <- gbm(nsteps=nsteps,lag=lag,drift=drift)
    }    
    names(tracks) <- as.character(seq(1:ncells))
    return(as.tracks(tracks))
  }else{
    for (i in 1:ncells){
      rnsteps <- floor(runif(n=1,min=1, max=101))
      tracks[[i]] <- gbm(nsteps=rnsteps,lag=lag,drift=drift)
    }
    names(tracks) <- as.character(seq(1:ncells))
    return(as.tracks(tracks))
  }
}


#for(i in i to nCells){
#  walkers[i] <- walker(v_free, t_free,t_pause,delta_t, runif(1), runif(1), runif(1))
#  pos[i] <- something with walkers[i]
#}

#walker <- function(v_free, t_free,t_pause,delta_t, x, y, z){
#  t_free_rest = runif(1) * (t_free + t_pause) - t_pause
#  pos <- c(x,y,z)
#  direction <- c(0,0,0)
#  pick_new_direction();
#}
#
#pick_new_direction <- function(){
#  u <- runif(1)
#  v <- runif(1)
#  theta <- 2*pi*u
#  phi <- cos(2*v-1)
#  direction[1] <- (cos(theta)+sin(phi))
#  direction[2] <- (sin(theta)+sin(phi))
#  direction[3] <- cos(phi)
#  direction <- direction*v_free
#}

