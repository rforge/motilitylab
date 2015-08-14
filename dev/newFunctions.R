

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


#' Simulates a single cell track based on model written by J. Textor, M. Sinn, and R. Boer called the Beauchemin model. 
#' 
#' @param nsteps Specifies the number of steps for the model to take
#' @param deltat Change in time between each timepoint
#' @param ppersist Indicates how probable a change in direction is. ppersist = 1 indicates there is never a change in direction between steps and ppersist = 0 indicates there is always a change in direction between steps
#' @param pbias Strength of movement in a given direction. Equates to the radius of sphere where the new direction is sampled. A pbias = 1 indicates sampling from unit sphere. 
#' @param taxis_mode Specified mode of movement. 1 := orthotaxis, 2 := topotaxis, 3 := klinotaxis
#' @param tfree Time interval for how long the cell is able to move between steps
#' @param vfree Speed of cell 
#' @param tpause Time for how long the cell takes to adjust movement to new direction
#' 
#' @details The default parameters are those found by mean-square displacement plots created by Beauchemin et al.
#' @return A matrix of size Nx5 where N is the number of timepoints allowed based number of \code{nsteps} allowed. The 
#' first column indicates the time change from the point previous to the current point. The second column is an 
#' index of timepoints (there are multiple timepoints allowed in one step). The remaining three columns are the 
#' three dimensional coordinates of the cell at that timepoint. 
#' @examples 
#' ## Create track with model parameters and return matrix of positions
#' out <- beauWalker(nsteps=20,ppersist = 0.3,taxis_mode = 1)
#' ## Plot X-Y projection
#' plot(PosY~PosX,data=out,type='b')
beauWalker <- function(nsteps = 10,deltat=1,ppersist=0.5,pbias=.9,taxis_mode=1,tfree=2,vfree=18.8,tpause=0.5){  
  #' Parameter checks
  if(ppersist < 0 || ppersist > 1){
    stop("ppersist must be a value between 0 and 1")
  }
  if(sum(taxis_mode == c(1,2,3))==0){
    stop("taxis_mode can either be 1, 2, or 3. See documentation for details on which model you would like to use.")
  }
  if(pbias <= 0){
    stop("pbias must be a value greater than 0")
  }
  if(nsteps <= 0){
    stop("nsteps must be a value larger than 0")
  }
  if(deltat <= 0){
    stop("deltat must be a value larger than 0")
  }
  if(tpause <= 0){
    stop("tpause must be a value larger than 0")
  }
  if(tfree <= 0){
    stop("tfree must be a value larger than 0")
  }
  if(vfree <= 0){
    stop("vfree must be a value larger than 0")
  }
  
  #' Helper function of beauWalker()
  #' This function samples uniformly from unit sphere
  sphere_sampler <- function(){
    x <- rnorm(3)
    z <- sqrt(sum(x^2))
    out <- c(x[1]/z,x[2]/z,x[3]/z)
    return(out)
    
  }
  
  #' Helper function of beauWalker()
  #' This function samples uniformly from unit circle
  circle_sampler <- function(){
    x <- rnorm(2)
    z <- sqrt(sum(x^2))
    out <- c(x[1]/z,x[2]/z)
    return(out)
  }
  
  #' Helper function of beauWalker()
  #' This function samples from the trianglular distribution
  triangular_sampler <- function(){
    x <- seq(from=-1,to=1,by=0.0001)
    y <- (x+1)/2
    z <- sample(x,1,prob=y/10000)
    return(z)
  }
  
  #' Helper function of beauWalker() and is not directly called by user
  pick_new_direction <- function(pbias,taxis_mode,tfree,vfree){
    d <- sphere_sampler()
    
    if(taxis_mode == 2){
      if(runif(1) < pbias){
        # topotaxis
        d[3] <- triangular_sampler()
        circ <- circle_sampler()
        circle_scale <- sqrt(1-d[3]^2)
        d[1] <- circ[1]*circle_scale
        d[2] <- circ[2]*circle_scale
      }
      d <- c(d*vfree,0)
    }
    if(taxis_mode == 3){
      if(pbias > 0){
        # klinotaxis
        tfreerest = tfreerest + sum(tfree * pbias * d);
      }
      d <- c(d*vfree,tfreerest)  
    }
    else{
      # orthotaxis
      d <- c(d*vfree*(1+pbias*d),0)
    }
    return(d)
  }
  
  #' Helper function of beauWalker() and is not directly called by user
  pausedStep <- function(deltat){
    if(-tfreerest+deltat < tpause ){
      tfreerest = tfreerest-deltat
      newpos <- c(i*deltat*60, i, pos[length(pos[,1]),3:5])
      return(c(tfreerest,1,newpos))
    }
    else{
      deltatrest= deltat - tpause - tfreerest
      tfreerest = tfree
      paused = FALSE
      
      if (runif(1) < 1-ppersist){
        d <- pick_new_direction(pbias=pbias,taxis_mode=taxis_mode,tfree=tfree,vfree=vfree)
        if(d[4]==0){
          d <- d[1:3]
        }
        else{
          tfreerest <- d[4]
          d <- d[1:3]
        }
      }
      return(c(deltatrest,2,d,tfreerest))
    }
  }
  
  #' Helper function of beauWalker() and is not directly called by user
  movingStep <- function(deltat){
    if( deltat < tfreerest ){
      newpos <- c(i*deltat*60,i,pos[length(pos[,1]),3:5] + d*deltat)
      tfreerest = tfreerest-deltat
      return(c(tfreerest,1,newpos))
    }
    else{
      deltatrest = deltat - tfreerest
      newpos <- c(i*deltat*60,i,pos[length(pos[,1]),3:5] + d*tfreerest)
      paused = TRUE
      tfreerest = 0.0 
      return(c(deltatrest,2,newpos))
    }
  }
  
  #' Initializing parameters 
  pos <- matrix(c(deltat,0,0,0,0),1,5)
  tfreerest = tfree
  paused = tfreerest < 0.0
  deltatresti <- deltatrest <- deltat - tpause - tfreerest
  d <- pick_new_direction(pbias=pbias,taxis_mode=taxis_mode,tfree=tfree,vfree=vfree)
  if(d[4]==0){
    d <- d[1:3]
  }
  else{
    tfreerest <- d[4]
    d <- d[1:3]
  }
  
  i <- 1
  j <- 0
  
  
  #' This is the primary loop that creates the 3d track
  repeat{
    if(paused==TRUE){
      catch <- pausedStep(deltat)
      if(catch[2]==1){
        pos <- rbind(pos,catch[3:7])
        i <- i+1
        tfreerest <- catch[1]
      }
      else{
        deltatrest <- catch[1]
        d <- catch[3:5]
        tfreerest<-catch[6]
        paused = FALSE
      }
      if( deltatrest > 0 ){
        catch <- movingStep(deltatrest)
        if(catch[2]==1){
          tfreerest <- catch[1]
          pos <- rbind(pos,catch[3:7])
          i <- i+1
        }
        else{
          deltatrest <- catch[1]
          pos <- rbind(pos,catch[3:7])
          i <- i+1
        }
      }
    }
    else{
      catch <- movingStep(deltat)
      if(catch[2]==1){
        tfreerest <- catch[1]
        pos <- rbind(pos,catch[3:7])
        i <- i+1
      }
      else{
        deltatrest <- catch[1]
        pos <- rbind(pos,catch[3:7])
        i <- i+1
        paused <- TRUE
      }
      if( deltatrest > 0 ){
        catch <- pausedStep(deltatrest)
        if(catch[2]==1){
          tfreerest <- catch[1]
          pos <- rbind(pos,catch[3:7])
          i <- i+1
        }
        else{
          deltatrest <- catch[1]
        }
      }
    }
    j <- j+1
    #or should we use i here instead? it's more intuitive and predictable if we use i... 
    if(j == nsteps){
      break
    }
  }
  colnames(pos) <- c("Delta_t","TimePoint","PosX","PosY","PosZ")
  return(pos)
}