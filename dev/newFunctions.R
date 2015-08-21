#' Helper function of beauWalker() and is not directly called by user
#' This function samples uniformly from unit sphere
.sphereSampler <- function(){
  x <- rnorm(3)
  z <- sqrt(sum(x^2))
  out <- x/z
  return(out)
}

#' Helper function of beauWalker() and is not directly called by user
#' This function samples uniformly from unit circle
.circleSampler <- function(){
  x <- rnorm(2)
  z <- sqrt(sum(x^2))
  out <- x/z
  return(out)
}

#' Helper function of beauWalker() and is not directly called by user
#' This function samples from the trianglular distribution
.triangularSampler <- function(){
  x <- runif(1)
  out <- -1+sqrt(4*x)
  return(out)
}


#' Helper function of beauWalker() and is not directly called by user
#' returns a direction for the cell to travel based on model parameters specified in beauWalker()
.pickNewDirection <- function(p.bias,taxis.mode,t.free,v.free,t.free.rest){
  d <- .sphereSampler()
  
  if(taxis.mode == 2){
    if(runif(1) < p.bias){
      # topotaxis
      d[3] <- .triangularSampler()
      circ <- .circleSampler()
      circle.scale <- sqrt(1-d[3]^2)
      d[1] <- circ[1]*circle.scale
      d[2] <- circ[2]*circle.scale
    }
    d <- c(d*v.free,0)
  }
  else if(taxis.mode == 3){
    if(p.bias > 0){
      # klinotaxis
      t.free.rest = t.free.rest + sum(t.free * p.bias * d);
    }
    d <- c(d*v.free,t.free.rest)  
  }
  else{
    # orthotaxis
    d <- c(d*v.free*(1+p.bias*d),0)
  }
  return(d)
}


#' Helper function of beauWalker() and is not directly called by user
#' Returns the next step in the beauWalker() if the cell is paused
.pausedStep <- function(d,delta.t,t.free.rest,t.pause,p.bias,taxis.mode,t.free,v.free,pos,i,p.persist){
  if(-t.free.rest+delta.t < t.pause ){
    t.free.rest = t.free.rest-delta.t
    new.pos <- c(i, pos)
    return(c(t.free.rest,1,new.pos))
  }
  else{
    delta.t.rest= delta.t - t.pause - t.free.rest
    t.free.rest = t.free
    
    if (runif(1) < 1-p.persist){
      d <- .pickNewDirection(t.free.rest=t.free.rest,p.bias=p.bias,taxis.mode=taxis.mode,t.free=t.free,v.free=v.free)
      if(d[4]==0){
        d <- d[1:3]
      }
      else{
        t.free.rest <- d[4]
        d <- d[1:3]
      }
    }
    return(c(delta.t.rest,2,d,t.free.rest))
  }
}


#' Helper function of beauWalker() and is not directly called by user
#' Returns the next step in the beauWalker() if the cell is moving
.movingStep <- function(delta.t,t.free.rest,d,pos,i){
  if( delta.t < t.free.rest ){
    new.pos <- c(i,pos + d*delta.t)
    t.free.rest = t.free.rest-delta.t
    return(c(t.free.rest,1,new.pos))
  }
  else{
    delta.t.rest = delta.t - t.free.rest
    new.pos <- c(i,pos + d*t.free.rest)
    return(c(delta.t.rest,2,new.pos))
  }
}

#' Simulates a single cell track based on updated Beauchemin model written by J. Textor, M. Sinn, and R. Boer. 
#' 
#' @param n.steps Specifies the number of steps for the model to take.
#' @param delta.t Change in time between each timepoint.
#' @param p.persist Indicates how probable a change in direction is. p.persist = 1 indicates there is never a 
#' change in direction between steps and p.persist = 0 indicates there is always a change in direction between steps.
#' @param p.bias Strength of movement in a given direction. Equates to the radius of sphere by which a new direction
#' is sampled. A p.bias = 1 indicates sampling from unit sphere. 
#' @param taxis.mode Specified mode of movement. 1 := orthotaxis, 2 := topotaxis, 3 := klinotaxis.
#' @param t.free Time interval for how long the cell is allowed to move between steps.
#' @param v.free Speed of cell .
#' @param t.pause Time for how long the cell takes to adjust movement to new direction.
#' 
#' @details The default parameters are those found by mean-square displacement plots created by Beauchemin et al.
#' @return A matrix of size \code{n.steps}x5. Each row represents the three dimensional coordinates of the cell
#' at the given timepoint.
#' @examples 
#' ## Create track with model parameters and return matrix of positions
#' out <- beauWalker(n.steps=20,p.persist = 0.3,taxis.mode = 1)
#' ## Plot X-Y projection
#' plot(y~x,data=out,type='b',main="Planar projection over z-axis")
beauWalker <- function(n.steps = 25,delta.t=1,p.persist=0.5,p.bias=0.9,taxis.mode=1,t.free=2,v.free=18.8,t.pause=0.5){  
  #' Parameter checks
  if(p.persist < 0 || p.persist > 1){
    stop("p.persist must be a value between 0 and 1")
  }
  if(sum(taxis.mode == c(1,2,3))==0){
    stop("taxis.mode can either be 1, 2, or 3. See documentation for details on which model you would like to use.")
  }
  if(p.bias < 0){
    stop("p.bias must be greater than or equal to 0")
  }
  if(n.steps <= 0){
    stop("n.steps must be larger than 0")
  }
  if(delta.t <= 0){
    stop("delta.t must be larger than 0")
  }
  if(t.pause < 0){
    stop("t.pause must be 0 or larger")
  }
  if(t.free <= 0){
    stop("t.free must be larger than 0")
  }
  if(v.free <= 0){
    stop("v.free must be larger than 0")
  }
  
  #' Initializing parameters 
  pos <- matrix(c(0,0,0,0),1,4)
  t.free.rest = t.free
  paused = t.free.rest < 0.0
  delta.t.rest <- delta.t - t.pause - t.free.rest
  d <- .pickNewDirection(t.free.rest=t.free.rest,p.bias=p.bias,taxis.mode=taxis.mode,t.free=t.free,v.free=v.free)
  if(d[4]==0){
    d <- d[1:3]
  }
  else{
    t.free.rest <- d[4]
    d <- d[1:3]
  }
  i <- 1
  
  #' This is the primary loop that creates the 3d track
  repeat{
    if(paused==TRUE){
      catch <- .pausedStep(d=d,p.persist=p.persist,i=i,delta.t=delta.t,t.free.rest=t.free.rest,t.pause=t.pause,p.bias=p.bias,taxis.mode=taxis.mode,t.free=t.free,v.free=v.free,pos=pos[length(pos[,1]),2:4])
      if(catch[2]==1){
        pos <- rbind(pos,catch[3:6])
        i <- i+1
        t.free.rest <- catch[1]
      }
      else{
        delta.t.rest <- catch[1]
        d <- catch[3:5]
        t.free.rest<-catch[6]
        paused = FALSE
      }
      if( delta.t.rest > 0 ){
        catch <- .movingStep(i=i,delta.t=delta.t.rest,pos=pos[length(pos[,1]),2:4],t.free.rest=t.free.rest,d=d)
        if(catch[2]==1){
          t.free.rest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
        else{
          t.free.rest <- 0.0
          delta.t.rest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
      }
    }
    else{
      catch <- .movingStep(i=i,delta.t=delta.t,pos=pos[length(pos[,1]),2:4],t.free.rest=t.free.rest,d=d)
      if(catch[2]==1){
        t.free.rest <- catch[1]
        pos <- rbind(pos,catch[3:6])
        i <- i+1
      }
      else{
        t.free.rest <- 0.0
        delta.t.rest <- catch[1]
        pos <- rbind(pos,catch[3:6])
        i <- i+1
        paused <- TRUE
      }
      if( delta.t.rest > 0 ){
        catch <- .pausedStep(d=d,p.persist=p.persist,i=i,delta.t=delta.t.rest,t.free.rest=t.free.rest,t.pause=t.pause,p.bias=p.bias,taxis.mode=taxis.mode,t.free=t.free,v.free=v.free,pos=pos[length(pos[,1]),2:4])
        if(catch[2]==1){
          t.free.rest <- catch[1]
          pos <- rbind(pos,catch[3:6])
          i <- i+1
        }
        else{
          delta.t.rest <- catch[1]
        }
      }
    }
    if(i >=n.steps+1){
      break
    }
  }
  colnames(pos) <- c("t","x","y","z")
  pos[,2:4]<- round(pos[,2:4],digits = 3)
  return(pos)
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