#' @title slidingWindowForStochasticMaps
#'
#' @description This function performs a sliding window analysis to obtain the number of species present in each region within regular time intervals based on the window size. This allows for comparison across biogeographic stochastic maps.
#' The default window size is 0.1 time units, but users should select values that make sense for their particular datasets.
#' @param tr Phylogenetic tree used in BSM analysis.
#' @param LTST.dt Lineages through space and time data table.
#' @param window Sets the window size. default = 0.1.
#' @return Returns an LTST data table with standardised event times.
#' @author Alex Skeels
#' @importFrom phytools nodeHeights
#' @export

slidingWindowForStochasticMaps <- function(LTST.dt, window=0.1){

  # standardise node hieghts so root = 1
  
  crown <- LTST.dt$time[1]
  # windows at set up from the root to the tips at intervals of the window size
  windows <- seq(from=0, to=crown, by=window)

  # create data frame to store stanardised time slice values
  window.df <- data.frame(matrix(0,ncol=length(LTST.dt), nrow=length(windows)))
  colnames(window.df) <- colnames(LTST.dt)
  window.df$time <- windows

  # iterate through windows
  for(i in 1:length(window.df[,1])){

    # find all event times that fall within the curent window
    current.window.rows <- which(LTST.dt$time >= window.df$time[i] & LTST.dt$time < window.df$time[i+1])

    # if no events fall in window then the current window inherets the values of the previous window
    if (length(current.window.rows) == 0 ) {
      
      if(i == 1){window.df[i, 2:length(LTST.dt)] <- tail(LTST.dt, 1)[2:length(LTST.dt)] }
      if(i > 1){window.df[i, 2:length(LTST.dt)] <- window.df[i -1, 2:length(LTST.dt)]}
    } else {
    # otherwise sum across columns for all rows in the current window
      window.df[i, 2:length(LTST.dt)] <- apply(LTST.dt[current.window.rows,2:length(LTST.dt)], MARGIN=2, max)
    }
  }
  window.df
}
