set_up_spatial_decontamination <- function(P4, well_dists){
  is_contaminant <- row.names(P4) %>% str_detect('extractNTC')
  SCR <- list()
  COVERAGE <- 10000
  background_contam_init <- rescale(1 + P4[which(is_contaminant), ] %>% apply(MARGIN = 2, min) ) %>% t()
  
  for(i in 1:nrow(P4)){
    neighbors <- which( (well_dists[i,] ==1 )&(well_dists[i,] > 0) )
    
    sources <- rbind(P4[neighbors,], background_contam_init)
    rare_sink <- 1    
    if(is_contaminant[i]==0){
      SCR[[i]] <- list( 
        neighbors_idx=neighbors,
        observed=P4[i,]
      )
    }else{
      q <- nrow(sources)
      SCR[[i]]=list(
        neighbors_idx=neighbors,
        observed=P4[i,]
      )
    }
  }
  return(SCR)
}


get_positions_df <- function(P4){
  wells <- row.names(P4) %>%
  sapply(function(x){
    q <- str_split(x, '[.]')[[1]]
    return( q[length(q)] %>% substr(1,1))
  }) %>%
  sapply( function(a) which(LETTERS == a) )


indices <- row.names(P4) %>%
  sapply(function(x){
    q <- str_split(x, '[.]')[[1]]
    return( q[length(q)] %>% substr(2,4))
  }) %>% as.integer()



positions <- data.frame(wells, indices) %>%
  mutate(indices=as.integer(indices)) %>% 
  as.matrix()  

row.names(positions) <- row.names(P4)

positions <- data.frame(positions)
positions$is_negative_control <- positions %>% row.names() %>% sapply(function(x) str_split(x, '[.]')[[1]][3]) == 'extractNTC'
return(positions)
}
