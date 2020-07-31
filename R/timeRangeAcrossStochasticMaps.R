#' @title timeRangeAcrossStochasticMaps
#'
#' @description Summarises the upper (97.5%) and lower (2.5%) quantiles of species diversity count at each time interval across biogeographic stochastic maps (BSMs).
#' The summarised data can be returned in 'melted' data frame with five columns: time, region, lowerQ, uperrQ, mid.
#' Alternatively the summarised data can be returned in a data frame with a column for time, and a further two columns for each geographic region summarising the upper and lower quantile of diversity at each time point.
#' @param tr Phylogenetic tree used in BSM analysis.
#' @param window.LTST Lineages through space and time data table with standardised event times (from slidingWindowForStochasticMaps).
#' @param melt Whether to melt the results (useful for plotting with ggplot2)
#' @return Diversity data for each geographic region summarised as quantiles in a data frame.
#' @author Alex Skeels
#' @import reshape2
#' @export


timeRangeAcrossStochasticMaps <- function(tr,window.LTST, melt=T){

require(reshape2)

  # count the number of unique geographic areas
  region.number <- length(window.LTST[[1]][, 2:length(window.LTST[[1]])])

  if(melt==F){
    # iterates over each column in the LTST data table that contains a geographic region.
    # for each geographic region, we create a data frame with each row representing a time point and each column is a different LTST data table from the window.LTST list
    # We can then apply the function quantile over each row (MARGIN=1) to get the upper a lower quantiles for each time point.
    # We store the results in a seperate column for upper and lower quantiles for each region in the time.range.df
    time.range.df <- data.frame(matrix(0, ncol= (1+ 2*region.number), nrow=length(window.LTST[[1]][,1])))
    names(time.range.df) <- c("time", sapply(colnames(window.LTST[[1]][,2:(1+region.number)]), FUN=function(x){paste(x, "max", collapse="", sep="_")} ), sapply(colnames(window.LTST[[1]][,2:(1+region.number)]), FUN=function(x){paste(x, "min", collapse="", sep="_")} ))
    for(i in 2:length(window.LTST[[1]])){
      region.name <- names(window.LTST[[1]])[i]
      region.diversity <- do.call("cbind", lapply(window.LTST, FUN=function(x){x[,i]}))
      region.diversity.max <- apply(region.diversity, MARGIN=1,quantile, probs=0.975 )
      time.range.df[, which(colnames(time.range.df) %in% paste(region.name, "upperQ", sep="_", collapse=""))] <- region.diversity.max
      region.diversity.min <- apply(region.diversity, MARGIN=1,quantile, probs=0.025 )
      time.range.df[, which(colnames(time.range.df) %in% paste(region.name, "lowerQ", sep="_", collapse=""))] <- region.diversity.min
    }
  } else {

    # if melt == TRUE
    # iterates over each column in the LTST data table that contains a geographic region.
    # for each geographic region, we create a data frame with each row representing a time point and each column is a different LTST data table from the window.LTST list
    # We can then apply the function quantile over each row (MARGIN=1) to get the upper a lower quantiles for each time point.
    # however differently from melt==FALSE, we store the results in a data frame with 4 columns and assign rows in the data frame to different regions using rows.to.occupy as an index.
    # this way data for all regions are stored in the same columns which makes it easy to visualise the results by plotting with ggplot2.
    # at the end we add an additional column "mid" as the mid point between the upper and lower quantile which may also be useful for plotting.
    time.range.df <- data.frame(matrix(0, ncol= 4, nrow=length(window.LTST[[1]][,1])*region.number))
    names(time.range.df) <- c("time", "region", "lowerQ", "upperQ")


    for(i in 2:length(window.LTST[[1]])){
      rows.to.occupy <- c((((i-1) * length(window.LTST[[1]][,1])) - length(window.LTST[[1]][,1])+1): ((i-1) * length(window.LTST[[1]][,1])))

      region.name <- names(window.LTST[[1]])[i]
      region.diversity <- do.call("cbind", lapply(window.LTST, FUN=function(x){x[,i]}))
      region.diversity.max <- apply(region.diversity, MARGIN=1,quantile, probs=0.975 )
      region.diversity.min <- apply(region.diversity, MARGIN=1,quantile, probs=0.025 )


      time.range.df[rows.to.occupy, "region"] <- region.name
      time.range.df[rows.to.occupy, "upperQ"] <- round(region.diversity.max)
      time.range.df[rows.to.occupy, "lowerQ"] <- round(region.diversity.min)
    }

    time.range.df$mid <- (time.range.df$lowerQ + time.range.df$upperQ) /2
  }

  time.range.df$time <- window.LTST[[1]]$time
  time.range.df
}
