#' @title getLTSTDataTable
#'
#' @description Creates a lineages through space and time data table by amalgamating the running diversity count across state combination identifiers into their component geographic areas.
#' @param ETT The events timing table (from getEventsTiming).
#' @param rangeStates The lookup table for geographic range states (from getRangeState).
#' @return An LTST data table as class 'data frame'
#' @author Alex Skeels
#' @import stringr
#' @export

getLTSTDataTable <- function(ETT, rangeStates){
  # requires stringr to identify unique character states
  require(stringr)
  unique.biogeo <- unique(unlist(stringr::str_extract_all(rangeStates$range, boundary("character"))))

  # creates new ETT except this time with just the unique geographic states as columns rather than state combination identifiers
  ETT2 <- as.data.frame(matrix(ncol = (length(unique.biogeo)+1), nrow=length(ETT[,1])))
  colnames(ETT2) <- c("time", unique.biogeo)
  ETT_tmp <- ETT
  colnames(ETT_tmp) <- rangeStates$range[match(colnames(ETT), rangeStates$states)]
  # loop to populate geographic states at each time step by macthing the original ETT to the lookup table and tallying the absolute diversity for each geographic state
  for(j in 1:length(unique.biogeo)){
    ETT2[, which(colnames(ETT2) %in% unique.biogeo[j])] <- rowSums(ETT_tmp[4:length(ETT_tmp)][grepl(unique.biogeo[j], colnames(ETT_tmp)[4:length(ETT_tmp)])])
  }
  ETT2$time <- ETT$time
  ETT2
}
