#' @title getEventsTiming
#'
#' @description Extracts the events timing information from the anagenetic and Cladogenetic events tables from a BGB analysis and stores it in the events timing table.
#' The events timing table contains the time of each event,
#' the node at which the event is located on the phylogeny (or beneath which for an anagenetic event along a branch),
#' whether the event was cladogenetic or anagenetic, the number of species extant at that time,
#' and a running total of the number of species present in each unique state combination.
#'
#' @param tree The phylogeny used in BSM analysis.
#' @param CET The cladogenetic events table output from BSM analysis.
#' @param AET The anagenetic events table output from BSM analysis.
#' @return Returns the events timing table as a data frame.
#' @author Alex Skeels
#' @export

getEventTiming <- function(tree, CET, AET){

  CET_k <- CET
  AET_k <-AET
  # get a vector of all named character states
  states <- unique(c(CET_k$sampled_states_AT_nodes, AET_k$current_rangenum_1based))
  # setUpETT creates a event timing table (ETT)
  ETT <- setUpETT(states, CET_k)
  # set the current event time to the root node and index which row this is in the CET
  event.time <- ETT[1,1]
  index <- which(CET_k$node.type == "root")

  # iterate through each event along the branches of the phlogeny forward through time (not including the tips)
  for (j in 1:((length(CET_k[,1])+length(AET_k[,1])))){


    # is the current event cladogenetic? - returns logical (TRUE for cladogenetic or FALSE for anagenetic)
    clado <- isEventClado(event.time, CET_k, AET_k)

    # finds the row number from the CET or AET corresponding to the current event
    index.node <- getTableIndexFromEventTime(CET_k, AET_k, ET =event.time, event.type=clado, ETT)
    index <- index.node[1]

    # specifies in the ETT whether the event is cladogenetic or anagenetic
    ETT[length(ETT[,1]), "event"] <- clado

    # if the current event is along a terminal branch (tip) in the phylogeny - select the next event time and move ahead if cladogenetic
    # otherwise identify which states have changed in the current time step
    # character changes can be anagenetic (lineage changes states with no net gain in diversity) or cladogenetic (lineage divides with a net increase in diversity)
    # whichStateChangesClado / whichStateChangesAna make changes to the event timing table directly
    if(clado == "clado"){ if(CET_k$node.type[index] == "tip"){event.time <- nextEventTime(CET_k, AET_k, event.time);next} else {ETT <- whichStateChangesClado(CET_k, ETT, event.time, index)}}
    if(clado == "ana"){ if(AET_k$node.type[index] == "tip"){ETT <- whichStateChangesAna(AET_k, ETT, event.time, index)}}

    # place extra data in ETT for current state - species diversity and node number
    ETT[length(ETT[,1]), "node"] <- index.node[2]
    ETT[length(ETT[,1]), "extant_sp"] <- sum(ETT[length(ETT[,1]), 5:length(ETT)])

    # select next chronological event time
    event.time <- nextEventTime(CET_k, AET_k, event.time)

    # add new row to ETT which will be altered in next iteration
    ETT <- rbind(ETT, ETT[length(ETT[,1]),], make.row.names=F)
    ETT[length(ETT[,1]), "time"] <- event.time
  }
  ETT
}

#' @title setUpETT
#'
#' @description creates the initial (unfilled) events timing table from a vector of character states
#' @param states vector of named character states
#' @param CET the cladogenetic events table output from BGB analysis
#' @return an empty ETT data frame
#' @author Alex Skeels
#' @export

setUpETT <- function(states, CET){

  # create data table with 4 columns for time, node, event, and diversity (extant_sp), as well as an additional column for each unique state
  state.no <- length(states)
  ETT <- data.frame(time=NA, "node"=0, "event" = NA, "extant_sp"=0)
  for( j in 1:state.no){ETT <- cbind(ETT, states[j]) }
  names(ETT)[5:length(ETT)] <- as.character(states)
  ETT[1,5:length(ETT)]<- rep(0, times=state.no)
  # populate ETT with values from the root node of the phylogeny
  root.node <- CET[CET$node.type == "root", "node"]
  root.value <- CET[CET$node.type == "root","sampled_states_AT_nodes"]
  ETT[1, 1] <- CET[CET$node.type == "root","time_bp"]
  ETT[1,"node"] <- root.node

  # return ETT
  ETT
}

#' @title isEventClado
#'
#' @description asks whether the current event is cladogenetic or anagenetic
#' @param event.time timing of the current event
#' @param CET the cladogenetic events table output from BGB analysis
#' @param AET the anagenetic events table output from BGB analysis
#' @return logical (TRUE for cladogenetic or FALSE for anagenetic)
#' @author Alex Skeels
#' @export


isEventClado <- function(event.time, CET, AET){

  if(any(CET$time_bp == event.time)){ return("clado") }

  if(any(AET$time_bp == event.time)){ return("ana") }
}


#' @title getTableIndexFromEventTime
#'
#' @description for either the CET or AET finds the index (row number) and node number for the current event using the event time as a reference
#' @param CET the cladogenetic events table output from BGB analysis
#' @param AET the anagenetic events table output from BGB analysis
#' @param event.type clado or ana
#' @param ETT events timing table
#' @return logical (TRUE for cladogenetic or FALSE for anagenetic)
#' @author Alex Skeels
#' @export

getTableIndexFromEventTime <- function(CET, AET, ET, event.type, ETT){

  # for either the CET or AET finds the index (row number) and node number for the current event
  if(event.type == "clado"){

    if(length(which(CET$time_bp %in% ET))==1){
      index <- which(CET$time_bp %in% ET)
      node <- CET$node[index]

    } else {
      index <- which(CET$time_bp %in% ET & !CET$node %in% ETT$node)[1]
      node <- CET$node[index]
    }

  }

  if(event.type == "ana"){

    if(length(which(AET$time_bp == ET))==1) {
      index <- which(AET$time_bp == ET)
      node <- AET$node[index]

    } else {

      index <- which(AET$time_bp == ET & !AET$node %in% ETT$node)[1]
      node <- AET$node[index]
    }
  }

  return(c(index, node))
}

#' @title whichStateChangesClado
#'
#' @description amodifies ETT by accounting for state changes occur at nodes in the phylogeny - cladogenetic events.
#' A state can increase in diversity by 1 if new daughter species also occurs in the same state.
#' A different state from the parent can increase by 1 if new daughter species occupies a new state.
#' A state can decrease by 1 if both daughter species occupy different states from the parent - in which case those new states will each +1
#' @param x the CET
#' @param ETT events timing table
#' @param ET current event time
#' @param index index of CET
#' @return updated ETT with changes made to the geographic states based on which event is occuring
#' @author Alex Skeels
#' @export


whichStateChangesClado <- function(x, ETT, ET, index){

  # $sampled_states_AT_brbots is the state number at the bottom of the branch below the node
  # $samp_RIGHT_dcorner is the state number of the right side brnach from the current node
  # $samp_LEFT_dcorner is the state number of the left side brnach from the current node

  #check that there are descendents from cladoegentic event
  if(is.na(x$samp_LEFT_dcorner[index]) & is.na(x$samp_RIGHT_dcorner[index])) {return(ETT)}

  # sampled_states_AT_brbots is the state number at the bottom of the branch below the node
  # will be NA for the root
  # if it is the root then add 1 to each state above the root (first two branches of the phylogeny)
  if(is.na(x$sampled_states_AT_brbots[index])) {
    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))] <- state.value + 1

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))] <- state.value + 1

    return(ETT)
  }

  # for all other branches:
  # 1 if descendants are both the same as ancestor (here we add 1 to the state of the ancestor)
  if(x$sampled_states_AT_brbots[index] == x$samp_LEFT_dcorner[index] & x$sampled_states_AT_brbots[index] == x$samp_RIGHT_dcorner[index]){

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$sampled_states_AT_brbots[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$sampled_states_AT_brbots[index]))] <- state.value + 1
    return(ETT)
  }

  #2 if the descendents are different from ancestor (here we minus 1 from ancestor character state and add 1 to  each new character state)
  if(x$sampled_states_AT_brbots[index] != x$samp_LEFT_dcorner[index] & x$sampled_states_AT_brbots[index] != x$samp_RIGHT_dcorner[index]){

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$sampled_states_AT_brbots[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$sampled_states_AT_brbots[index]))] <- state.value - 1

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))] <- state.value + 1

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))] <- state.value + 1
    return(ETT)
  }

  #3  if one descendant is different from ancestor (here we add 1 to new character state)
  if(x$sampled_states_AT_brbots[index] != x$samp_LEFT_dcorner[index] & x$sampled_states_AT_brbots[index] == x$samp_RIGHT_dcorner[index]){

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_LEFT_dcorner[index]))] <- state.value + 1
    return(ETT)
  }

  if(x$sampled_states_AT_brbots[index] == x$samp_LEFT_dcorner[index] & x$sampled_states_AT_brbots[index] != x$samp_RIGHT_dcorner[index]){

    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$samp_RIGHT_dcorner[index]))] <- state.value + 1

    return(ETT)

  }

}


#' @title whichStateChangesAna
#'
#' @description modifies ETT by accounting for state changes along branches in the phylogeny (anagenetic changes).
#' The values of the originally occupied state will decrease by 1 and values of the new state will increase by 1.
#' @param x the AET
#' @param ETT events timing table
#' @param ET current event time
#' @param index index of AET
#' @return updated ETT with changes made to the geographic states based on which event is occuring
#' @author Alex Skeels
#' @export

whichStateChangesAna <- function(x, ETT, ET, index){

  # make sure the anagenetic change is to a new state
  if(x$current_rangenum_1based[index] != x$new_rangenum_1based[index]) {

    # deduct one from original state
    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$current_rangenum_1based[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$current_rangenum_1based[index]))] <- state.value - 1
    # add one to new state
    state.value <- ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$new_rangenum_1based[index]))]
    ETT[which(ETT[,1] == ET), which(names(ETT) == as.character(x$new_rangenum_1based[index]))] <- state.value + 1

    return(ETT)
  }
}

#' @title nextEventTime
#'
#' @description this function takes the current event time and finds which event time folows it chronological by indexing the AET and CET
#' @param CET cladogenetic event table
#' @param AET anagenetic event table
#' @param ET current event time
#' @return returns the next event time in chronological order
#' @author Alex Skeels
#' @export

nextEventTime <- function(CET, AET, ET) {

  event.timing <- c(CET$time_bp, AET$time_bp)

  order.event <- event.timing[order(event.timing)]

  nextEventTime <- order.event[which(order.event == ET)- 1]

  if(order.event[which(order.event == ET)] == head(order.event)[1]){ nextEventTime=0}
  return(nextEventTime)
}


