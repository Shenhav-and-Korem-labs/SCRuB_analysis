install.packages("mvtnorm", repos='http://cran.us.r-project.org')



library(tidyverse)
library(magrittr)
library(vegan)
library(philentropy)

source('Simulation_Functions/Plate_simulations.R')
source('Simulation_Functions/generate_simulations.R')

source('SCRuB/lsq_initializations.R')
source('SCRuB/main_functions.R')
source('SCRuB/spatial_functions.R')

set.seed(1)

run_just_spatial_SCRUB <- function(samples, controls, well_dist, dist_threshold=1){
  combined <- rbind(samples, controls)
  is_contaminant <- c( rep(F, nrow(samples)),  rep(T, nrow(controls)))
  spatial_SCRUB_out <- spatial_SCRUB( combined, is_contaminant,well_dist, dist_threshold=dist_threshold)
  return(  list(spatial_SCRUB_out=spatial_SCRUB_out
  ) ) 
}

## load SparseDOSSA2 functiond
setwd('Fig4/SparseDOSSA2/R')
list.files()[1:13] %>% sapply(source)
setwd('../../..')

# load the data simulators
load('../data/Fig4/Fig_4_simulators.Rdata')

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

store_results_with_sample_type_and_jsd <- function(new_results, 
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
                                           sample_type, 
                                           control_type,
                                           idx){
  
  
  input_sink_jsd <- get_rowwise_jsd(samples, 
                                       get_rescaled_mat(observed_samples))
  
  spatial_SCRUB_jsd <- get_rowwise_jsd(samples, 
                                     get_rescaled_mat( out$spatial_SCRUB_out$decontaminated_samples ) )
  
  jsds <- c(input_sink_jsd, spatial_SCRUB_jsd)
  names <-  c(rep('Input', nrow(observed_samples)), 
              rep('Spatial SCRUB', nrow(observed_samples)))
  
  result <- data.frame(names, jsds) %>% 
    mutate( n_taxa=n_taxa, 
            n_controls=n_controls, 
            contam_species=sum(true_contamination>0),
            contam_level=contam_level, 
            well_to_well_level=well_level, 
            plate_format=plate_format, 
            well_species_limit=well_species_limit,
            sample_type=sample_type,
            control_type=control_type,
            round=idx) 
  if(is.data.frame(new_results)==F){ return(result)
  }else{return(rbind(new_results, result))}
}


run_fig4_simulations <- function(n_taxa=100, 
                                 n_samples=95, 
                                 control_noise='med', 
                                 sample_noise='med', 
                                 contam_in_sample_noise='med',
                                 n_experiments_per_set = 100,
                                 contam_range=c(.15), 
                                 well_to_well_range=c(.1), 
                                 plate_format_range=c('A', 'B', 'C', 'D', 'E', 'F'), 
                                 n_control_range=c(2,4,8), 
                                 exp_type='Fig_4'
                                 ){

total_experiments <- n_experiments_per_set * 3 * length(plate_format_range)*length(n_control_range) * length(contam_range) * length(well_to_well_range)
full_fig4_results <- NA
idx <- 1

for( sample_type in c('high', 'med', 'low')){
  if(sample_type=='high')sample_generator <- fitted_sparsedossa_PNP_samples
  if(sample_type=='med')sample_generator <- fitted_hands
  if(sample_type=='low')sample_generator <- fitted_vag_samples
  
  for( control_type in c('low')){
    if(control_type=='low')control_generator <- fitted_cancer_extractions     
    for(contam_level in contam_range){
      for(n_contam_species in c(NA)){
        for(experiment in 1:n_experiments_per_set){
          
          n_controls <- 8
          
          samples <- SparseDOSSA2(template = sample_generator, 
                                  n_sample=n_samples, 
                                  n_feature = n_taxa, 
                                  median_read_depth = 10000, 
                                  verbose = F)$simulated_data %>% t() %>% get_rescaled_mat()

          row.names(samples) <- paste0('Sample_', 1:n_samples)

          true_contamination <- (SparseDOSSA2(template = control_generator, 
                                              n_sample=2, 
                                              n_feature = n_taxa, 
                                              median_read_depth = 10000, 
                                              verbose = F)$simulated_data %>% t())[1, ] %>% rescale()
          
          
          if(is.na(n_contam_species)==F){
            true_contamination[order(true_contamination, decreasing = T)[1+n_contam_species:length(true_contamination)] ] <- 0
            true_contamination <- rescale(true_contamination)
          }
          

          # draw noisy / contaminated samples, 8 noisy controls
          observed_samples <- sapply(1:n_samples, function(x){
            draw_contaminated_sample(samples[x,],
                                     true_contamination, 
                                     samples,
                                     contam_level, 
                                     n_taxa, 
                                     c(),
                                     sample_noise=sample_noise, 
                                     sample_well_level=0,
                                     contam_in_sample_noise=contam_in_sample_noise
                                     ) %>% return() 
                                     }  
          ) %>% t()
          
          row.names(observed_samples) <- row.names(samples)
          
          
          controls_pre_well_to_well <- sapply(1:n_controls, function(x){
            draw_control( true_contamination, 
                          samples, 
                          x, 
                          c(), 
                          0,
                          control_noise, 
                          n_taxa
            ) %>% return()
          } ) %>% t()
          
          row.names(controls_pre_well_to_well) <- paste0('Control_', 1:8)
          
          
          for( well_level in well_to_well_range){
            for(plate_format in plate_format_range){
              for(n_controls in n_control_range){
              n_smps <- 96-n_controls
              well_positions <- Simulate_Plate_Positions(n_smps, n_controls, n_wells=8, plate_format=plate_format)
              well_dists <- dist(well_positions) %>% as.matrix()
              well_species_limit <- NA
              if(n_controls>1){ tmp <- controls_pre_well_to_well[1:n_controls,]
              }else{
                tmp <- controls_pre_well_to_well[1, ] %>% t()
                row.names(tmp) <- c('Control_1')
                }
              
              observed_controls <- get_final_conts(tmp,
                                                   samples[1:n_smps,],
                                                   well_dists, 
                                                   well_level, 
                                                   well_species_limit)
                if(n_controls>1){ out <- run_just_spatial_SCRUB(observed_samples[1:n_smps,], observed_controls[1:n_controls,], well_dists, dist_threshold=1.5)
                }else{
                  tmp <- observed_controls[1, ] %>% t()
                  row.names(tmp) <- c('Control_1')
                  out <- run_just_spatial_SCRUB(observed_samples[1:n_smps,], tmp, well_dists, dist_threshold=1.5)
                }
               full_fig4_results <- store_results_with_sample_type_and_jsd(
                                                              full_fig4_results, 
                                                               out,
                                                                samples[1:n_smps,], 
                                                  observed_samples[1:n_smps,], 
                                                              n_controls,
                                                          well_species_limit, 
                                                             plate_format, 
                                                         well_level, 
                                                             contam_level, 
                                                             true_contamination,
                                                             n_taxa,
                                                                sample_type,
                                                               control_type,
                                                           idx)

                print(paste('Running', exp_type, 'iteration', 
                  idx, 'of', total_experiments))
                idx <- idx+1
              }
            }
          }
        }
      }
    }
    
  }
}

return(full_fig4_results)
}

full_f4_results <- run_fig4_simulations(n_experiments_per_set=5) # set to 100 to repeat original analysis
dir.create('../results/data/Fig_4')
full_f4_results %>% write.csv('../results/data/Fig_4/Fig_4_Simulations.csv')             
               

# simulations varying contamination, well contamination
full_s7_results <- run_fig4_simulations(n_experiments_per_set=1, 
                                       contam_range=c(.05, .25, .5), 
                                       well_to_well_range=c(.05, .1, .15, .2, .25, .3, .35, .4, .45, .5), 
                                       plate_format_range=c('A'), 
                                       n_control_range=c(2,4,8), 
                                       exp_type='Fig_S7') # set to 100 to repeat original analysis

full_s7_results %>% write.csv('../results/data/Fig_4/Fig_S7_Simulations.csv')
               
               
               
               
               


