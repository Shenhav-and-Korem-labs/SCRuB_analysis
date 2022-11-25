install.packages("Pareto", repos='http://cran.us.r-project.org')
library(Pareto)


rDirichlet <- function(n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

make_noise <- function(a, N){
  num_noises <-  3 
  noise_idx <- sample(N, size = num_noises, replace = FALSE)
  noise = rmultinom(1, a, rPareto(N, 1, 1) ) #general_noise
  noise[noise_idx] <- rmultinom(1, a, rescale( rPareto(num_noises, 1, 1) ) )   #specific noise
  return(noise)
}

all_noise <-  function(samp, level=25){
  level <- rpois(1, level) + round( sum(samp) / 50 )
  N  <- length(samp)
  rand_noise <- make_noise(25, N) 
  resample <- function(cont) rmultinom( 1, level,  rescale(cont))
  x = ( resample(samp) + rand_noise ) * ( 2 * rbinom(N, 1, .5) - 1 )
  #ensure we don't observe negative reads
  x[ (samp + x) < 0 ] = -samp[ (samp+x) < 0]
  return( x  + samp)
}

make_rand_mixing <- function(K, mu){
  props <- rnorm(1, mean=mu, sd=.1)
  props[props>.95] <- .95
  props[props<.1] <- .1
  return( generate_mixing(K, props))
}

generate_sample_dists <- function(n_taxa, 
                                  sparsity){
  
  if(sparsity=='high') n_draws=5
  if(sparsity=='med') n_draws=50
  if(sparsity=='low') n_draws=200
  
  return( make_noise(n_draws, n_taxa) %>% rescale )
}


simulate_samples <- function(n_taxa, 
                             n_samples, 
                             sample_similarity){
  
  if(sample_similarity=='high') r <- .3
  if(sample_similarity=='med') r <- .75
  if(sample_similarity=='low') r <- 1
  
  q <- rPareto(n_taxa, .1, r) 
  q[q> 15] <- 15
  q <- q %>% rescale()
  w <- rDirichlet(n_samples, q)
  return(w)
}


draw_contaminated_sample <- function(s_dist, 
                                     c_dist, 
                                     samples,
                                     contam_level, 
                                     n_taxa,
                                     well_contaminants_vec, 
                                     sample_noise='med', 
                                     contam_in_sample_noise='med',
                                     sample_well_level=0.05,
                                     COVERAGE=10000){
  if(sample_noise=='high') s_noise = .05 *(  make_noise(100, n_taxa) %>% rescale() )
  if(sample_noise=='med') s_noise = .01 * ( make_noise(10, n_taxa ) %>% rescale() ) 
  if(sample_noise=='low') s_noise = .0005 * ( make_noise(10, n_taxa) %>% rescale() )
  if(sample_noise=='none') s_noise = 0
  
  if(contam_in_sample_noise=='high') c_noise = .05 * ( make_noise(100, n_taxa) %>% rescale() )
  if(contam_in_sample_noise=='med') c_noise = .01 * ( make_noise(10, n_taxa) %>% rescale() ) 
  if(contam_in_sample_noise=='low') c_noise=  .0005 * ( make_noise(10, n_taxa)%>% rescale() )
  if(contam_in_sample_noise=='none') c_noise= 0
  
  
  if(contam_in_sample_noise!='none') c_noise[c_noise>.15] <- 0
  if(sample_noise!='none') s_noise[s_noise>.15] <- 0

  neighbors <- names(well_contaminants_vec[ which( well_contaminants_vec == 1 ) ] )
  well_amt = rnorm(1, sample_well_level, .02)
  if(well_amt<0) well_amt <- 0
  if(well_amt>1) well_amt <- 1
  if(sample_well_level==0) well_amt <- 0
  
  if(length(neighbors)>1){ 
    mix <- make_rand_mixing( length(neighbors)-1, .8 ) %>% sample
    well_contam_profile <- ( t(samples[neighbors,  ]) * mix ) %>% rowSums()
  }else if(length(neighbors)==1){well_contam_profile <- samples[neighbors,]
  }else{well_contam_profile <- rep(0, length(s_dist) ) }
  
  tmp <- ( well_contam_profile*well_amt ) +
    ( 1-well_amt )*(  (1-contam_level)*rescale(s_dist+s_noise[1:length(s_dist)]) + 
                        contam_level*rescale(c_dist +c_noise[1:length(c_dist)]) )
                      
  tmp[tmp<0] <- 0

  rmultinom( 1, COVERAGE, tmp
               ) %>%
    return()
}

simulate_well_dist <- function(n_controls, n_samples){
  well_dists <- matrix(10, nrow = n_samples + n_controls, ncol = n_samples + n_controls)
  for(i in 1:n_controls){
    neighbors <- sample(n_samples, 4, replace=F)
    well_dists[n_samples + i, neighbors] <- 1
  }
  return(well_dists)
}


draw_control <- function(true_contaminant,
                         samples,
                         control_idx,
                         well_contaminants,
                         well_to_well_level,
                         control_noise='med',
                         n_taxa,
                         COVERAGE=10000
){
  
  if(control_noise=='high') c_noise <-  .1 * ( make_noise(500, n_taxa) %>% rescale() )
  if(control_noise=='med') c_noise <-  .01 * ( make_noise(10, n_taxa)%>% rescale() ) 
  if(control_noise=='low') c_noise <-  .0005 * ( make_noise(10, n_taxa)%>% rescale() )
  if(control_noise=='none') c_noise <- 0
  if(control_noise!='none') c_noise[c_noise>.2] <- 0
  
  well_contaminants_vec <- well_contaminants[paste0('Control_', control_idx),]
  neighbors <- c( names(well_contaminants_vec[ which( well_contaminants_vec == 1 ) ] ) )
  well_amt <- rnorm(1, well_to_well_level, .04)
  if(well_amt<0)well_amt <- 0
  if(well_amt>1)well_amt <- 1
  if(well_to_well_level==0) well_amt <- 0
  if(length(neighbors)>1){ 
    mix <- make_rand_mixing( length(neighbors)-1, .8 ) %>% sample
    well_contam_profile <- ( t(samples[neighbors,  ]) * mix ) %>% rowSums()
  }else if(length(neighbors)==1){well_contam_profile <- samples[neighbors,]
  }else{well_contam_profile <- rep(0, length(true_contaminant) ) }
  
  
  c_noise <-  c_noise[1:length(true_contaminant)]
  c_noise[is.na(c_noise)] <- 0

  rmultinom( 1, COVERAGE, (well_amt)*well_contam_profile +  (1-well_amt)*( rescale(true_contaminant + c_noise ) )) %>%
    return()
}

set_up_run <- function( n_taxa = 100,
                        n_samples = 50,
                        n_controls = 5,
                        sample_similarity = 'med',
                        control_sparsity = 'med',
                        sample_noise='med', 
                        contam_in_sample_noise='med',
                        control_noise = 'med', 
                        contam_level = .2,
                        well_to_well_level = .1, 
                        use_real_data=F, 
                        sample_well_level=0.05,
                        plate_format='A'){
  
  if(use_real_data==F){
    samples <- SparseDOSSA2(template = fitted_sparsedossa_samples, 
                            n_sample=n_samples, 
                            n_feature = n_taxa, 
                            median_read_depth = 10000, 
                            verbose = F)$simulated_data %>% t() %>% get_rescaled_mat()
    
    
    true_contamination <- (SparseDOSSA2(template = fitted_sparsedossa_controls, 
                                        n_sample=2, 
                                        n_feature = n_taxa, 
                                        median_read_depth = 10000, 
                                        verbose = F)$simulated_data %>% t())[1, ] %>% rescale()
  }else{
    rows <- sample( nrow(samples), size = n_samples, replace = FALSE)
    samples <- samples[rows, ] %>% get_rescaled_mat()
    true_contamination <- controls[sample(nrow(controls), 1), ] %>% rescale()
    n_taxa <- ncol(samples)
  }
  
  row.names(samples) <- paste0('Sample_', 1:n_samples)
  
  plate_positions <- Simulate_Plate_Positions(n_samples, n_controls, plate_format = plate_format)
  well_dist <- plate_positions %>% dist(method='euclidean') %>% as.matrix()
  
  well_contaminants <- simulate_well_contaminants(well_dist)
  well_contaminants <- well_contaminants[,which( colnames(well_contaminants) %>% str_detect('Sample') )]
  
  observed_samples <- sapply(1:n_samples, function(x){
    draw_contaminated_sample(samples[x,],
                             true_contamination, 
                             samples,
                             contam_level, 
                             n_taxa, 
                             well_contaminants[x,],
                             sample_noise=sample_noise, 
                             sample_well_level=sample_well_level,
                             contam_in_sample_noise=contam_in_sample_noise) %>% return() }  
  ) %>% t()
  
  observed_controls <-sapply(1:n_controls, function(x){
    draw_control( true_contamination, 
                  samples, 
                  x, 
                  well_contaminants, 
                  well_to_well_level,
                  control_noise, 
                  n_taxa
    ) %>% return()
  } ) %>% t()
  
  row.names(observed_controls) <- paste0('Control_', 1:n_controls)
  row.names(observed_samples) <- paste0('Sample_', 1:n_samples)
  
  return(list(underlying_samples = samples,
              true_contamination=true_contamination, 
              well_dist = well_dist, 
              observed_samples=observed_samples, 
              observed_controls=observed_controls, 
              n_taxa=n_taxa, 
              n_samples=n_samples, 
              n_controls=n_controls, 
              sample_similarity=sample_similarity,
              control_sparsity=control_sparsity,
              sample_noise=sample_noise, 
              contam_in_sample_noise=contam_in_sample_noise, 
              control_noise=control_noise,
              contam_level=contam_level, 
              well_to_well_level=well_to_well_level, 
              plate_format=plate_format, 
              sample_well_level=sample_well_level, 
              plate_positions=plate_positions))
}



run_just_spatial_SCRUB <- function(samples, controls, well_dist, dist_threshold=1){
  combined <- rbind(samples, controls)
  is_contaminant <- c( rep(F, nrow(samples)),  rep(T, nrow(controls)))
  spatial_SCRUB_out <- spatial_SCRUB( combined, is_contaminant,well_dist, dist_threshold=dist_threshold)
  return(  list(spatial_SCRUB_out=spatial_SCRUB_out
  ) ) 
}



JSD <- function(x, test.na = TRUE, unit = "log2", est.prob = NULL){
  ## overwriting jsd function to mute messaging
  if (!is.matrix(x))
    stop("Please provide a matrix as input, e.g. with x <- rbind(vector1, vector2).", call. = FALSE)
  
  return( philentropy::distance(x = x, 
                   method   = "jensen-shannon", 
                   test.na  = test.na, 
                   unit     = unit,
                   est.prob = est.prob,
                   mute.message = TRUE) )
}

get_rowwise_jsd <- function(A, B){
  zz <- sapply(1:nrow(A), function(x) JSD( rbind( rescale(A[x, ]), rescale(B[x, ]) ) ) ) %>% unname() %>%
    return()
}



store_results <- function(new_results, 
                          out,
                          samples, 
                          observed_samples, 
                          n_controls, 
                          well_species_limit, 
                          plate_format, 
                          well_level, 
                          contam_level, 
                          true_contamination,
                          n_taxa,
                          idx){
  
  
  input_sink_jsd <- get_rowwise_jsd(samples, 
                                       get_rescaled_mat(observed_samples))
  
  spatial_SCRUB_jsd <- get_rowwise_jsd(samples, 
                                     get_rescaled_mat( out$spatial_SCRUB_out$decontaminated_samples ) )
  
  jsds <- c(input_sink_jsd, spatial_SCRUB_jsd)
  names <-  c(rep('Input', length(input_sink_jsd) ), 
              rep('Spatial SCRUB', length(input_sink_jsd) )
              )
              

 result <- data.frame(names, jsds) %>% 
    mutate( n_taxa=n_taxa, 
            n_controls=n_controls, 
            contam_species=sum(true_contamination>0),
            contam_level=contam_level, 
            well_to_well_level=well_level, 
            plate_format=plate_format, 
            well_species_limit=well_species_limit,
            round=idx) 
  if(is.data.frame(new_results)==F){ return(result)
  }else{return(rbind(new_results, result))}
}




simulate_well_contaminants <- function(well_dists){
  well_coefs<-1+well_dists*0#matrix of 1s with the sample shape as dist matrix
  
  # add a few samples that contaminate more often
  num_big_leakers <- rpois(1, 3.5)
  big_leakers <- sample(row.names(well_dists), num_big_leakers)
  
  for(leaker in big_leakers) well_coefs[, leaker] <- rpois(1,8) 
  
  q <- well_coefs*1/(2*well_dists**3 )
  q[q>.75] <- .75
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

create_well_contam_mixing <- function(samples, well_contaminants){
  n_well_contams <- sum(well_contaminants>0)
  if(n_well_contams==0)return(0*samples[1, ])
  if(n_well_contams==1)return(samples[which(well_contaminants==1),])
  
  mixing <- rDirichlet(1, rep(1/n_well_contams, n_well_contams))
  
  return( ( samples[which(well_contaminants==1), ] * as.vector(mixing) ) %>% colSums() )
}


combine_well_mixings <- function(cont, well_mix, well_mix_amt){
  COVERAGE <- 10000
  well_amt <- rnorm(1, well_to_well_level, .04)
  if(well_amt<0)well_amt <- 0
  if(well_amt>1)well_amt <- 1
  if(well_to_well_level==0) well_amt <- 0
  if(sum(well_mix)==0) well_amt <- 0
  
  rmultinom( 1, COVERAGE, (well_amt)*well_mix +  (1-well_amt)*(cont) ) %>%
    return()
}



get_final_conts <- function(controls_pre_well_to_well, samples, well_dists, well_amt, well_species_limit){
  
  row.names(controls_pre_well_to_well) <- row.names(well_dists)[which(str_detect(row.names(well_dists), 'Control'))]
  well_contaminants <- simulate_well_contaminants(well_dists)
  well_contaminants <- well_contaminants[which(str_detect(row.names(well_contaminants), 'Control')), 
                                         which(str_detect(row.names(well_contaminants), 'Sample'))]
  
  
  well_contam_mixings <- well_contaminants %>% 
    apply(MARGIN=1, function(x) return( create_well_contam_mixing(samples,x) ) ) %>% t()
  
  if(is.na(well_species_limit)==F){
    if(ncol(samples) > well_species_limit){
      mixings <- well_contam_mixings %>% colSums()
      well_contam_mixings[,order(mixings, decreasing = T)[1+well_species_limit:length(mixings)] ] <- 0
      well_contam_mixings <- get_rescaled_mat(well_contam_mixings)
    }}
  well_contam_mixings[is.na(well_contam_mixings)] <- 0
  
  scaled_observed_conts <- controls_pre_well_to_well %>% get_rescaled_mat()
  
  final_observed_conts <- scaled_observed_conts*0
  
  for(cont in row.names(final_observed_conts)){
    final_observed_conts[cont,] <- combine_well_mixings(scaled_observed_conts[cont,],well_contam_mixings[cont,], well_amt )
  }
  
  return( final_observed_conts )
}


combine_well_mixings <- function(cont, well_mix, well_mix_amt){
  COVERAGE <- 10000
  well_amt <- rnorm(1, well_mix_amt, .04)
  if(well_amt<0)well_amt <- 0
  if(well_amt>1)well_amt <- 1
  if(well_mix_amt==0) well_amt <- 0
  if(sum(well_mix)==0) well_amt <- 0
  rmultinom( 1, COVERAGE, (well_amt)*well_mix + (1-well_amt)*(cont) ) %>%
    return()
}


