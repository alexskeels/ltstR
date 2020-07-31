#' @title getLTSTFromStochasticMaps
#'
#' @description Applies getLTSTDataTable across two lists of cladogenetic and anagenetic events tables from a biegoegraphic stochastic mapping analysis (BSM).
#'This function iterates over each stochastic map in the list, producing a list of LTST data tables.
#' @param tr Phylogenetic tree used in BGB analysis.
#' @param CET_list List of cladogenetic events tables from a BSM analysis.
#' @param AET_list List of Anagenetic events tables from a BSM analysis.
#' @return A list of LTST data tables.
#' @author Alex Skeels
#' @export

getLTSTFromStochasticMaps <- function(tr, CET_list, AET_list, tipranges){

  if(class(CET_list) != "list"){print("CET_list and AET need to be lists from stochastic mapping in BioGeoBEARs")}

  # create object to store bLTT curve data in and iterate over each AET/CET in the lists
  LTST.dt.list <- vector("list", length(CET_list))

  for(i in 1:length(CET_list)){

    CTT <- getEventTiming(tr, CET_list[[i]], AET_list[[i]])
    rangeStates <- getRangeState(CET_list[[i]], AET_list[[i]], tipranges)
    LTST.dt <- getLTSTDataTable(CTT, rangeStates)

    LTST.dt.list[[i]] <- LTST.dt
  }

  return(LTST.dt.list)
}
