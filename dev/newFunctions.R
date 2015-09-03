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