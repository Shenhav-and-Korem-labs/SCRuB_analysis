
Simulate_Plate_Positions <- function(n_samples, n_controls, n_wells=8, plate_format='A'){
  n_cols <- ceiling( (n_samples+n_controls)/n_wells )
  
  structure_A_neg_positions <- c(1, 1) %>%
    rbind( c(n_wells, n_cols)) %>%
    rbind( c(1, n_cols)) %>%
    rbind( c(n_wells, 1)) %>%
    rbind( c(1, floor(n_cols/2))) %>%
    rbind( c(n_wells, ceiling(n_cols/2)) ) %>%
    rbind( c(floor(n_wells/2), 1)) %>%
    rbind( c(ceiling(n_wells/2), n_cols))
  
  structure_B_neg_positions <- c(ceiling(n_wells/3), floor(n_cols/5)) %>%
    rbind( c(ceiling(2*n_wells/3), ceiling(3*n_cols/5)) ) %>%
    rbind( c(ceiling(2*n_wells/3), floor(n_cols/5)) ) %>%
    rbind( c(ceiling(n_wells/3), ceiling(3*n_cols/5)) ) %>%
    rbind( c(ceiling(n_wells/3), ceiling(4.5*n_cols/5)) ) %>%
    rbind( c(ceiling(2*n_wells/3), ceiling(4.5*n_cols/5)) ) %>%
    rbind( c(ceiling(2*n_wells/3), ceiling(2*n_cols/5)) ) %>%
    rbind( c(ceiling(n_wells/3), ceiling(2*n_cols/5)) ) 
  
  if(plate_format=='A') neg_spots <- structure_A_neg_positions[1:n_controls,]
  if(plate_format=='B') neg_spots <- structure_B_neg_positions[1:n_controls,]
  if(plate_format=='C') neg_spots <- cbind(1:n_wells, floor(n_cols/2) )
  if(plate_format=='D') neg_spots <-  cbind(1:n_wells, rep(1,n_wells) )
  if(plate_format=='E'){
    neg_tmp <- 1:3 %>% merge(1:3) %>% as.matrix()
    neg_spots <- neg_tmp[neg_tmp %>% apply(MARGIN=1, function(x)sum(x)+max(x)) %>% order(decreasing = FALSE), ][1:8, ]
  }
  if(plate_format=='F'){ combos <- expand.grid(x = seq(1, n_wells, 1), y = seq(1, n_cols, 1)) %>%
    tbl_df()
                 neg_spots <- sample_n(combos, size = 8, replace = FALSE) %>% as.matrix()
  }
  
  positions_df <- 1:n_wells %>% merge(1:n_cols)
  
  if(n_controls==1) neg_spots <- neg_spots %>% t()
  for(j in 1:n_controls){
    row.names(positions_df)[ (positions_df$x==neg_spots[j, 1])&(
      positions_df$y==neg_spots[j, 2]) ] <- c(paste0('Control_', j))
  }
  
  tmp <- positions_df[str_detect(row.names(positions_df), 'Control')==FALSE, ]
  row.names(tmp) <- 1:nrow(tmp)
  tmp <- tmp[as.double( row.names(tmp) )<=n_samples,]
  row.names(tmp) <-  paste0('Sample_', row.names(tmp) )
  positions_df <- rbind(positions_df[str_detect(row.names(positions_df), 'Control'),], tmp)
  return( positions_df[ row.names(positions_df)[ order(row.names(positions_df)) ], ] )
}




simulate_well_contaminants <- function(well_dists){
  q <- 1/(2*well_dists**3 )
  q[ is.infinite(q) ] <- 0
  
  determine_well_contams <- function(a){
    w <- c()
    for(i in a) w <- c(w, rbinom(1, 1, i) )
    names(w) <- names(a)
    return(w)
  }
  
  q %>% apply(MARGIN=1, determine_well_contams) %>%
    return()
}
