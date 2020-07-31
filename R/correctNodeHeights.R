#' @title correctNodeHeights
#'
#' @description adds a small random number to each node height so that none are matching
#' @param CAET List or singular cladogenetic events or anagenetic events tables from a BSM analysis.
#' @param sd Standard deiation of normal deviate to add to node heights
#' @return CAET
#' @author Alex Skeels
#' @export


correctNodeHeights <- function(CET, AET, sd=0.01){

    if(class(CET) == "data.frame"){

      CET$time_bp <- CET$time_bp + rnorm(length(CET$time_bp), 0, sd)
      AET$time_bp <- AET$time_bp + rnorm(length(AET$time_bp), 0, sd)

      while(any(duplicated(c(CET$time_bp, AET$time_bp)))){
        CET$time_bp <- CET$time_bp + rnorm(length(CET$time_bp), 0, sd)
        AET$time_bp <- AET$time_bp + rnorm(length(AET$time_bp), 0, sd)
      }
    }
    if(class(CET) == "list"){


      CET <- lapply(CET, FUN=function(x){x$time_bp <- x$time_bp + rnorm(length(x$time_bp), 0, sd); return(x)})
      AET <- lapply(AET, FUN=function(x){x$time_bp <- x$time_bp + rnorm(length(x$time_bp), 0, sd); return(x)})

      for(i in 1:length(CET)){

        while(any(duplicated(c(CET[[i]]$time_bp, AET[[i]]$time_bp)))){
          CET[[i]]$time_bp <- CET[[i]]$time_bp + rnorm(length(CET[[i]]$time_bp), 0, sd)
          AET[[i]]$time_bp <- AET[[i]]$time_bp + rnorm(length(AET[[i]]$time_bp), 0, sd)
        }

      }
    }


  return(list(CET, AET))

}
