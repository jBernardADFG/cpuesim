#' Function to sample from the simulated cpue according to a predefined sampling design
#' 
#' @param cpue_sim A data.frame containing the simulated cpue values and coordinate locations
#' @param n_sets The total number of sets
#' @return Returns a list of data.frames where each contains the transect #, set within transect #, and cpue value
#' @export
two.stage.sample <- function(cpue_sim, n_sets){
  sample_list <- list()
  for (sim in 3:ncol(cpue_sim)){
    cpue_sim_i <- cbind(cpue_sim[,1],cpue_sim[,2],cpue_sim[,sim])
    possible_transects <- unique(cpue_sim_i$x)
    n_sets <- 0
    samples <- data.frame()
    while(T){
      tran <- sample(possible_transects, 1)
      possible_transects <- possible_transects[possible_transects!=tran]
      sets <- cpue_sim_i[cpue_sim_i$x==tran, ] 
      if(nrow(samples)+nrow(sets) < n_sets){
        samples <- rbind(samples, sets)
      }else{
        samples <- rbind(samples, sets[1:(n_sets-nrow(samples)),])
        sample_list[[sim-2]] <- samples
        break
      }
    }
  }
  return(sample_list)
}

