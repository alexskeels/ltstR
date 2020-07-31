#' @title correctNodeHeights
#'
#' @description adds a small random number to each node height so that none are matching
#' @param CAET List or singular cladogenetic events or anagenetic events tables from a BSM analysis.
#' @param sd Standard deiation of normal deviate to add to node heights
#' @return CAET
#' @author Alex Skeels
#' @export


correctNodeHeights <- function(CAET, sd=0.01){
  
    if(class(CAET) == "data.frame"){
      
      CAET$time_bp <- CAET$time_bp + rnorm(length(CAET$time_bp), 0, sd)
      while(!length(CAET$time_bp)== length(unique(CAET$time_bp))){
        CAET$time_bp <- CAET$time_bp + rnorm(length(CAET$time_bp), 0, sd)        
      }
    }
    if(class(CAET) == "list"){
      
      CAET <- lapply(CAET, FUN=function(x){x$time_bp <- x$time_bp + rnorm(length(x$time_bp), 0, sd); return(x)})
      while(!all(unlist(lapply(CAET, FUN=function(x){length(x$time_bp)== length(unique(x$time_bp))})))){
        CAET <- lapply(CAET, FUN=function(x){x$time_bp <- x$time_bp + rnorm(length(x$time_bp), 0, sd); return(x)})
        
      }
        
    }

  
  return(CAET)
  
}
