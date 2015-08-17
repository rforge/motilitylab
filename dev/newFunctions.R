

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


#' Helper function of beauWalker() and is not directly called by user
#' This function samples uniformly from unit sphere
.sphereSampler <- function(){
  x <- rnorm(3)
  z <- sqrt(sum(x^2))
  out <- c(x[1]/z,x[2]/z,x[3]/z)
  return(out)
  
}

#' Helper function of beauWalker() and is not directly called by user
#' This function samples uniformly from unit circle
.circleSampler <- function(){
  x <- rnorm(2)
  z <- sqrt(sum(x^2))
  out <- c(x[1]/z,x[2]/z)
  return(out)
}

#' Helper function of beauWalker() and is not directly called by user
#' This function samples from the trianglular distribution
.triangularSampler <- function(){
  x <- runif(1)
  y <- -1+sqrt(4*x)
  return(y)
}


#' Helper function of beauWalker() and is not directly called by user
#' returns a direction for the cell to travel based on model parameters specified in beauWalker()
.pickNewDirection <- function(pBias,taxisMode,tFree,vFree,tFreeRest){
  d <- .sphereSampler()
  
  if(taxisMode == 2){
    if(runif(1) < pBias){
      # topotaxis
      d[3] <- .triangularSampler()
      circ <- .circleSampler()
      circleScale <- sqrt(1-d[3]^2)
      d[1] <- circ[1]*circleScale
      d[2] <- circ[2]*circleScale
    }
    d <- c(d*vFree,0)
  }
  else if(taxisMode == 3){
    if(pBias > 0){
      # klinotaxis
      tFreeRest = tFreeRest + sum(tFree * pBias * d);
    }
    d <- c(d*vFree,tFreeRest)  
  }
  else{
    # orthotaxis
    d <- c(d*vFree*(1+pBias*d),0)
  }
  return(d)
}


#' Helper function of beauWalker() and is not directly called by user
#' Returns the next step in the beauWalker() if the cell is paused
.pausedStep <- function(d,deltaT,tFreeRest,tPause,pBias,taxisMode,tFree,vFree,pos,i,pPersist){
  if(-tFreeRest+deltaT < tPause ){
    tFreeRest = tFreeRest-deltaT
    newpos <- c(i, pos)
    return(c(tFreeRest,1,newpos))
  }
  else{
    deltaTrest= deltaT - tPause - tFreeRest
    tFreeRest = tFree
    
    if (runif(1) < 1-pPersist){
      d <- .pickNewDirection(tFreeRest=tFreeRest,pBias=pBias,taxisMode=taxisMode,tFree=tFree,vFree=vFree)
      if(d[4]==0){
        d <- d[1:3]
      }
      else{
        tFreeRest <- d[4]
        d <- d[1:3]
      }
    }
    return(c(deltaTrest,2,d,tFreeRest))
  }
}


#' Helper function of beauWalker() and is not directly called by user
#' Returns the next step in the beauWalker() if the cell is moving
.movingStep <- function(deltaT,tFreeRest,d,pos,i){
  if( deltaT < tFreeRest ){
    newpos <- c(i,pos + d*deltaT)
    tFreeRest = tFreeRest-deltaT
    return(c(tFreeRest,1,newpos))
  }
  else{
    deltaTrest = deltaT - tFreeRest
    newpos <- c(i,pos + d*tFreeRest)
    return(c(deltaTrest,2,newpos))
  }
}

#' TODO: For some combinations of paramters the function gets stuck?? Not sure if it's my computer or not... 
#' Simulates a single cell track based on model written by J. Textor, M. Sinn, and R. Boer called the Beauchemin model. 
#' 
#' @param nsteps Specifies the number of steps for the model to take
#' @param deltaT Change in time between each timepoint
#' @param pPersist Indicates how probable a change in direction is. pPersist = 1 indicates there is never a change in direction between steps and pPersist = 0 indicates there is always a change in direction between steps
#' @param pBias Strength of movement in a given direction. Equates to the radius of sphere where the new direction is sampled. A pBias = 1 indicates sampling from unit sphere. 
#' @param taxisMode Specified mode of movement. 1 := orthotaxis, 2 := topotaxis, 3 := klinotaxis
#' @param tFree Time interval for how long the cell is able to move between steps
#' @param vFree Speed of cell 
#' @param tPause Time for how long the cell takes to adjust movement to new direction
#' 
#' @details The default parameters are those found by mean-square displacement plots created by Beauchemin et al.
#' @return A matrix of size Nx5 where N is the number of timepoints allowed based number of \code{nsteps} allowed. The 
#' first column indicates the time change from the point previous to the current point. The second column is an 
#' index of timepoints (there are multiple timepoints allowed in one step). The remaining three columns are the 
#' three dimensional coordinates of the cell at that timepoint. 
#' @examples 
#' ## Create track with model parameters and return matrix of positions
#' out <- beauWalker(nsteps=20,pPersist = 0.3,taxisMode = 1)
#' ## Plot X-Y projection
#' plot(PosY~PosX,data=out,type='b')
beauWalker <- function(nsteps = 25,deltaT=1,pPersist=0.5,pBias=0.5,taxisMode=1,tFree=2,vFree=18.8,tPause=0.5){  
  #' Parameter checks
  if(pPersist < 0 || pPersist > 1){
    stop("pPersist must be a value between 0 and 1")
  }
  if(sum(taxisMode == c(1,2,3))==0){
    stop("taxisMode can either be 1, 2, or 3. See documentation for details on which model you would like to use.")
  }
  if(pBias < 0){
    stop("pBias must be a value greater than or equal to 0")
  }
  if(nsteps <= 0){
    stop("nsteps must be a value larger than 0")
  }
  if(deltaT <= 0){
    stop("deltaT must be a value larger than 0")
  }
  if(tPause <= 0){
    stop("tPause must be a value larger than 0")
  }
  if(tFree <= 0){
    stop("tFree must be a value larger than 0")
  }
  if(vFree <= 0){
    stop("vFree must be a value larger than 0")
  }
  
  #' Initializing parameters 
  pos <- matrix(c(0,0,0,0),1,4)
  tFreeRest = tFree
  paused = tFreeRest < 0.0
  deltaTrest <- deltaT - tPause - tFreeRest
  d <- .pickNewDirection(tFreeRest=tFreeRest,pBias=pBias,taxisMode=taxisMode,tFree=tFree,vFree=vFree)
  if(d[4]==0){
    d <- d[1:3]
  }
  else{
    tFreeRest <- d[4]
    d <- d[1:3]
  }
  i <- 1
  
  #' This is the primary loop that creates the 3d track
  repeat{
    if(paused==TRUE){
      catch <- .pausedStep(d=d,pPersist=pPersist,i=i,deltaT=deltaT,tFreeRest=tFreeRest,tPause=tPause,pBias=pBias,taxisMode=taxisMode,tFree=tFree,vFree=vFree,pos=pos[length(pos[,1]),2:4])
      if(catch[2]==1){
        pos <- rbind(pos,catch[3:6])
        i <- i+1
        tFreeRest <- catch[1]
      }
      else{
        deltaTrest <- catch[1]
        d <- catch[3:5]
        tFreeRest<-catch[6]
        paused = FALSE
      }
      if( deltaTrest > 0 ){
        catch <- .movingStep(i=i,deltaT=deltaTrest,pos=pos[length(pos[,1]),2:4],tFreeRest=tFreeRest,d=d)
        if(catch[2]==1){
          tFreeRest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
        else{
          tFreeRest <- 0.0
          deltaTrest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
      }
    }
    else{
      catch <- .movingStep(i=i,deltaT=deltaT,pos=pos[length(pos[,1]),2:4],tFreeRest=tFreeRest,d=d)
      if(catch[2]==1){
        tFreeRest <- catch[1]
        pos <- rbind(pos,catch[3:6])
        i <- i+1
      }
      else{
        tFreeRest <- 0.0
        deltaTrest <- catch[1]
        pos <- rbind(pos,catch[3:6])
        i <- i+1
        paused <- TRUE
      }
      if( deltaTrest > 0 ){
        catch <- .pausedStep(d=d,pPersist=pPersist,i=i,deltaT=deltaTrest,tFreeRest=tFreeRest,tPause=tPause,pBias=pBias,taxisMode=taxisMode,tFree=tFree,vFree=vFree,pos=pos[length(pos[,1]),2:4])
        if(catch[2]==1){
          tFreeRest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
        else{
          deltaTrest <- catch[1]
        }
      }
    }
    if(i == nsteps+1){
      break
    }
  }
  colnames(pos) <- c("t","x","y","z")
  return(pos)
}