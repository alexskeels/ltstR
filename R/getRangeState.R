#' @title getRangeState
#'
#' @description Generates lookup table to match geographic range states (as charater strings) to state combination identifiers (as numeric).
#' @param CET The cladogenetic events table output from BSM analysis.
#' @param AET The anagenetic events table output from BSM analysis.
#' @param tipranges The tipranges object (class 'tipranges') from BGB analysis containing a presence absence matrix with species names as rows and geographic states in columns.
#' @return A lookup table as class 'data frame'.
#' @author Alex Skeels
#' @import BioGeoBEARS
#' @export

getRangeState <- function(CET, AET, tipranges){

  # get vector naming all unique state combination identifiers
  states <- unique(c(CET$sampled_states_AT_nodes, AET$current_rangenum_1based))
  # create range.states lookup table
  range.states <- data.frame(range=NA, states=states)
  tipranges2 <- tipranges@df

  # replace 1's in tipranges with state names (to keep track of later)
  for(i in 1:length(tipranges2)){
    tipranges2[which(tipranges2[,i] == "1"),i] <- names(tipranges2)[i]
  }

  # add column to tipranges which combines names of all ranges each species occupies (might exist in more than one range)
  tipranges2$range.name <- NA
  for(i in 1:length(tipranges2[,1])){
    tipranges2$range.name[i] <- paste(tipranges2[i,which(tipranges2[i,] != 0)], collapse ="")
  }
  # find match between state combination identifier in CET_tips$sampled_states_AT_nodes and the actual names of the geographic ranges
  CET_tips <- CET[which(CET$node.type == "tip"),]
  for(i in 1:length(states)){
    if(states[i] %in% CET_tips$sampled_states_AT_nodes) {
      tip <- CET_tips$label[which(CET_tips$sampled_states_AT_nodes %in% states[i])][1]
      range.states$range[which(range.states$states == states[i])] <- tipranges2$range.name[which(rownames(tipranges2) == tip)]
    } else { next}
  }

  # some range names will be NA because the infromation on them is in the AET
  # so do the same but using the anagenetic events table (which contains the rest of the infromation)
  for(j in 1:length(states)){

    if(!is.na(range.states$range[which(range.states$states == states[j])])){next}

    if(states[j] %in% AET$current_rangenum_1based){
      range.states[j,1] <- AET$current_rangetxt[which(AET$current_rangenum_1based == states[j])[1]]
    } else {

      if(states[j] %in% AET$new_rangetxt){
        range.states[j,1] <- AET$new_rangetxt[which(AET$new_rangenum_1based == states[j])[1]]

      } else {

        if(states[j] %in% CET$sampled_states_AT_nodes){
          states.txt <- CET$clado_event_txt[which(CET$sampled_states_AT_nodes == states[j])]
          range.states[j,1] <- strsplit(unique(states.txt[which(states.txt != "")])[1], "-")[[1]][1]
        }
      }
    }
  }
  range.states
}
