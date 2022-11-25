
library(tidyverse)
library(magrittr)
# library(phyloseq)
library(decontam)
library(vegan)
library(philentropy)
library(microDecon)

source('SCRuB/lsq_initializations.R')
source('SCRuB/main_functions.R')
source('SCRuB/spatial_functions.R')

set.seed(1)

synth_data <- read.csv('../data/Fig1/synth_combination.csv', row.names = 1) %>% as.matrix()
colnames(synth_data) <- paste0('OTU_', 1:ncol(synth_data))

controls <- synth_data[row.names(synth_data) %>% str_detect('BLANK') %>% 
                         which(), ]
samples <- synth_data[( F == row.names(synth_data) %>% str_detect('BLANK') ) %>% 
                        which(), ]


source('Simulation_Functions/Plate_simulations.R')
source('Simulation_Functions/generate_simulations.R')


CLEAN_SAMPLES_DECONTAM <- function(samples, 
                                   controls, 
                                   threshold = 0.1
){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  
  ps <- samples %>% rbind(controls)
     
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isContaminant(ps, neg=neg, threshold = threshold,
                       detailed = T,
                       method='prevalence'
  )
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == T] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$contaminant),
              cont_idx = out, 
              estimated_sources = controls) )
}


CLEAN_SAMPLES_DECONTAM_low_biomass <- function(samples, 
                                               controls, 
                                               threshold = 0.5){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)

  ps <- samples %>% rbind(controls)
       
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isNotContaminant(ps, neg, threshold = threshold, 
                          detailed = T)
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == FALSE] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$not.contaminant),
              cont_idx = out, 
              estimated_sources = controls) )
}

CLEAN_SAMPLES_MICRODECON <- function(smps, cnts){
    tmp <- data.frame( rbind( colnames(smps), cnts, smps ) %>% t() )
    tmp[, 2:( ncol(tmp) ) ]  <- as.numeric( as.matrix( tmp[, 2:( ncol(tmp) ) ] ) )
    decontaminated <- decon(data = tmp, numb.blanks=nrow(cnts), numb.ind= c(nrow(smps)), taxa=F)
     
    md_out <- smps*0
    md_out[,decontaminated$decon.table$V1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t()
    
    return(list(decontaminated_samples=md_out))
}


get_samps_from_upload <- function(uploaded_samples, n_samples){
  n_provided <- nrow(uploaded_samples)
  return(uploaded_samples[sample(n_provided, n_samples, replace=n_samples>n_provided),
                          ])
}


REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X = samples
  to_remove <- which( (controls %>% colSums() )  > 0  )
  X[,to_remove] <- 0
  return(X)
}


run_SCRUB_DEC_and_RESTR<- function(samples, controls, well_dist, dist_threshold=1.5){
  combined <- rbind(samples, controls)
  is_contaminant <- c( rep(F, nrow(samples)),  rep(T, nrow(controls)))
  spatial_SCRUB_out <- spatial_SCRUB( combined, is_contaminant,well_dist,            
                                               dist_threshold=dist_threshold)
  
  SCRUB_out <-SCRUB(samples, controls)
  decont_out <- CLEAN_SAMPLES_DECONTAM(samples, controls)
  decont_low_out <- CLEAN_SAMPLES_DECONTAM_low_biomass(samples, controls)
  restr_out <- REMOVE_ALL_CONTROL_SPECS(samples, controls)
  
  
  return(  list(spatial_SCRUB_out=spatial_SCRUB_out, 
                SCRUB_out=SCRUB_out, 
                decont_out=decont_out, 
                decont_low_out=decont_low_out,
                microdec_out=microdec_out,
                restr_out=restr_out
  ) ) 
}

REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X = samples
  to_remove <- which( (controls %>% colSums() )  > 0  )
  X[,to_remove] <- 0
  return(X)
}


run_SCRUB_DEC_and_RESTR<- function(samples, controls, well_dist, dist_threshold=1.5){
  combined <- rbind(samples, controls)
  is_control <- c( rep(F, nrow(samples)),  rep(T, nrow(controls)))
  spatial_SCRUB_out <- spatial_SCRUB(combined, 
                                     is_control,
                                     well_dist,            
                                     dist_threshold=dist_threshold
                                     )

  SCRUB_out <- SCRUB(samples, controls)
  decont_out <- CLEAN_SAMPLES_DECONTAM(samples, controls)
  decont_low_out <- CLEAN_SAMPLES_DECONTAM_low_biomass(samples, controls)
  restr_out <- REMOVE_ALL_CONTROL_SPECS(samples, controls)
  microdec_out <- CLEAN_SAMPLES_MICRODECON(samples, controls)
  
  return(  list(spatial_SCRUB_out=spatial_SCRUB_out, 
                SCRUB_out=SCRUB_out, 
                decont_out=decont_out, 
                decont_low_out=decont_low_out,
                restr_out=restr_out,
                microdec_out=microdec_out
  ) ) 
}

store_f1_results <- function(new_results, 
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
  
  SCRUB_jsd <- get_rowwise_jsd(samples, 
                                get_rescaled_mat( out$SCRUB_out$decontaminated_samples ) )

  dec_jsd <- get_rowwise_jsd(samples, 
                              get_rescaled_mat( out$decont_out$decontaminated_samples ) )
  
  
  microdec_jsd <- get_rowwise_jsd(samples, 
                                  get_rescaled_mat( out$microdec_out$decontaminated_samples ) )

  q <- get_rescaled_mat( out$decont_low_out$decontaminated_samples )
  q[is.na(q)] <- 0
  dec_low_jsd <- get_rowwise_jsd(samples, 
                                  q)

  restrictive_jsd <- get_rowwise_jsd(samples, 
                                      get_rescaled_mat( out$restr_out) )
  
  jsds <- c(input_sink_jsd, spatial_SCRUB_jsd, SCRUB_jsd, 
             dec_jsd, dec_low_jsd, restrictive_jsd, microdec_jsd)
  names <-  c(rep('Input', length(input_sink_jsd)), 
              rep('Spatial SCRUB', length(input_sink_jsd)), 
              rep('SCRUB', length(input_sink_jsd)), 
              rep('Decontam', length(input_sink_jsd)), 
              rep('Decontam (Low Biomass)', length(input_sink_jsd)), 
              rep('Restrictive', length(input_sink_jsd)),
              rep('microDecon', length(input_sink_jsd))
                  )

  result <- data.frame(names, jsds) %>% 
    mutate( n_taxa=n_taxa, 
            n_controls=n_controls, 
            contam_species=sum(true_contamination>0),
            contam_level=contam_level, 
            well_to_well_level=well_level, 
            plate_format=plate_format, 
            round=idx) 
  if(is.data.frame(new_results)==F){ return(result)
  }else{return(rbind(new_results, result))}
}


finalized_simulations <- function(full_samps, 
                                  full_conts, 
                                  use_uploaded_samples=T,
                                  n_experiments=1){

  total_experiments <- n_experiments * 3 * 4 * 4

  n_taxa <- ncol(full_samps)
  n_samples <- 88
  control_noise <- 'high'
  sample_noise <- 'high'
  contam_in_sample_noise <- 'high'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  for(contam_level in c(.05, .25, .50) ){
    for(experiment in 1:n_experiments){
      
      n_controls <- 8
      
      
      samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
      row.names(samples) <- paste0('Sample_', 1:n_samples)
      
      true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
      n_taxa <- ncol(samples)
      
      
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
                                 contam_in_sample_noise=contam_in_sample_noise) %>% return() }  
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
      
      
      for(plate_format in c('F')){ 
        well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
        well_dists <- dist(well_positions) %>% as.matrix()
        for(well_level in c(0, .05, .25, .50)){
          well_species_limit <- NA
          observed_controls <- get_final_conts(controls_pre_well_to_well,
                                               samples,
                                               well_dists, 
                                               well_level, 
                                               well_species_limit)
          colnames(observed_samples) <- colnames(samples)
          colnames(observed_controls) <- colnames(samples)
          for(n_controls in c(8, 4, 2, 1)){ # test across 1:8 controls
            if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, observed_controls[1:n_controls,], well_dists, dist_threshold=1.5)
            }else{
              tmp <- observed_controls[1, ] %>% t()
              row.names(tmp) <- c('Control_1')
              out <- run_SCRUB_DEC_and_RESTR(observed_samples, tmp, well_dists, dist_threshold=1.5)
            }
            # store_finalized_results_using_jsd
            print(paste('Running Fig_1 iteration', 
                  idx, 'of', total_experiments))
            sitch_results <- store_f1_results(sitch_results, 
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
                                          idx)
            
            idx <- idx+1
          }
        }
      }
      print(idx)
    }
  }
  return(sitch_results)
}



results_of_final_sims <- finalized_simulations(samples, controls, n_experiments = 10) #10 for original analysis
# dir.create('../results/data/')
# dir.create('../results/data/Fig_1')

results_of_final_sims %>%
  write.csv('../results/data/Fig_1/F1_Simulation_Results.csv')



vary_well_level_sims <- function(full_samps, 
                                  full_conts, 
                                  use_uploaded_samples=T,
                                  n_experiments=1){
  
  total_experiments <- n_experiments * 2 * 3 * 4
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  control_noise <- 'med'
  sample_noise <- 'med'
  contam_in_sample_noise <- 'med'
  
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  for(contam_level in c(.05, .25) ){
    for(experiment in 1:n_experiments){
      
      n_controls <- 8
      
      
      samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
      row.names(samples) <- paste0('Sample_', 1:n_samples)
      
      true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
      n_taxa <- ncol(samples)
      
      
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
                                 contam_in_sample_noise=contam_in_sample_noise) %>% return() }  
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
      
      
      for(plate_format in c('F')){ 
        well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
        well_dists <- dist(well_positions) %>% as.matrix()
        for(well_level in c(.05, .25, .50)){
          well_species_limit <- NA
          observed_controls <- get_final_conts(controls_pre_well_to_well,
                                               samples,
                                               well_dists, 
                                               well_level, 
                                               well_species_limit)
          colnames(observed_samples) <- colnames(samples)
          colnames(observed_controls) <- colnames(samples)
          for(n_controls in c(8, 4, 2, 1)){ # test across 1:8 controls
            if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, observed_controls[1:n_controls,], well_dists, dist_threshold=1.5)
            }else{
              tmp <- observed_controls[1, ] %>% t()
              row.names(tmp) <- c('Control_1')
              out <- run_SCRUB_DEC_and_RESTR(observed_samples, tmp, well_dists, dist_threshold=1.5)
            }
            # store_finalized_results_using_jsd
            print(paste('Running Fig_1 iteration', 
                        idx, 'of', total_experiments))
            sitch_results <- store_f1_results(sitch_results, 
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
                                              idx)
            
            idx <- idx+1
          }
        }
      }
      print(idx)
    }
  }
  return(sitch_results)
}


set.seed(1)
results_of_welll_sims <- vary_well_level_sims(samples, controls, n_experiments = 10) 

results_of_welll_sims %>% write.csv('../results/data/Fig_1/F1_well_leakage_Results.csv')



run_s2_simulations <- function(samples, full_conts, n_experiments=50){
  
  run <- 1
  
  jsds <- c()
  names <- c()
  rounds <- c()
  smp_tps <- c()
  cont_tps <- c()
  input_jsd <- c()
  n_contam_taxa <- c()
  sample_cont_overlap <- c()
  n_samples <- 88
  n_controls <- 8
  
  removed_taxa=list()
  removed_taxa[['SCRUB']] <- list()
  removed_taxa[['dec']] <- list()
  removed_taxa[['dec (lb)']] <- list()
  
  cont_tp <- 1
  
  samp_tp='dorm'
  cont_tp <- 'extraction'
  
  for(test in 1:n_experiments){
    
    
    sample_simulation <- samples[sample(nrow(samples), 88 ), ]
    row.names(sample_simulation) <- paste0('Sample_', 1:n_samples)
    
    true_contamination <- full_conts[sample(nrow(full_conts), 1), ] %>% rescale
    n_taxa <- ncol(sample_simulation)
    
    control_simulation <- sapply(1:n_controls, function(x){
      draw_control( true_contamination, 
                    sample_simulation, 
                    x, 
                    c(), 
                    0,
                    'med', 
                    n_taxa
      ) %>% return()
    } ) %>% t()
    
    row.names(control_simulation) <- paste0('Control_', 1:8)
    
    row.names(sample_simulation) <- paste0('Sample_', 1:nrow(sample_simulation))
    
    well_dists <- Simulate_Plate_Positions(88, 8, 8, 'A') %>% dist('euclidean') %>% as.matrix()
    tmp_scr <- spatial_SCRUB( rbind(sample_simulation, control_simulation), 
                              is_control = c( rep(F, nrow(sample_simulation)), 
                                              rep(T, nrow(control_simulation))), 
                              well_dists = well_dists)
    
    
    dec <- CLEAN_SAMPLES_DECONTAM(sample_simulation, control_simulation )
    
    decontammed <- CLEAN_SAMPLES_DECONTAM_low_biomass(sample_simulation, control_simulation )
    
    microdeconned <-  CLEAN_SAMPLES_MICRODECON(sample_simulation, control_simulation)
    
    
    tmp_scr$decontaminated_samples[ rowSums(tmp_scr$decontaminated_samples)==0 , ] = 1
    decontammed$decontaminated_samples[ rowSums(decontammed$decontaminated_samples)==0 , ] = 1
    mean_controls <- control_simulation %>% colSums
    
    input_jsd <- c( input_jsd, sample_simulation %>% apply(MARGIN=1, function(x) get_rowwise_jsd(x %>% t(), mean_controls %>% t())) )
    
    jsds <- c( jsds, get_rowwise_jsd(sample_simulation, tmp_scr$decontaminated_samples), get_rowwise_jsd(sample_simulation, decontammed$decontaminated_samples), 
               get_rowwise_jsd(sample_simulation, dec$decontaminated_samples), get_rowwise_jsd(sample_simulation, microdeconned$decontaminated_samples) )
    
    
    names <- c(names, rep('SCRUB', 88), rep('Decontam (LB)', 88), rep('Decontam', 88), rep('microDecon', 88))
    rounds <- c(rounds, rep(run, 88*4))
    smp_tps <- c(smp_tps, rep(samp_tp, 88*4))
    cont_tps <- c(cont_tps, rep(cont_tp, 88*4))
    n_contam_taxa <- c(n_contam_taxa, sum( mean_controls > 0))
    sample_cont_overlap <- c(sample_cont_overlap, sample_simulation %>% apply(MARGIN=1, function(x) sum( (x %>% rescale)*(mean_controls %>% rescale) ) ) )
    
    run <- run+1
    print(run)
    
    
    # keep track of taxa added/removed
    started_with_idx <- ( sample_simulation %>% colSums() > 0 ) %>% unname() %>% which
    removed_taxa[['SCRUB']][[test]] <- ( (tmp_scr$decontaminated_samples %>% colSums() ==0 ) %>% unname() )[started_with_idx]
    removed_taxa[['dec (lb)']][[test]] <- ( (decontammed$decontaminated_samples %>% colSums() == 0 )%>% unname() )[started_with_idx]
    removed_taxa[['dec']][[test]] <- ( (dec$decontaminated_samples %>% colSums() == 0 )%>% unname() )[started_with_idx]
    removed_taxa[['microDecon']][[test]] <- ( (microdeconned$decontaminated_samples %>% colSums() == 0 )%>% unname() )[started_with_idx]
    
    
    
  }
  ## write jsd dataframe
  data.frame(jsds, names, rounds, smp_tps, cont_tps, input_jsd, n_contam_taxa, sample_cont_overlap)  %>%
    write.csv('../results/data/Fig_1/no_contamination_simulation.csv')
  
  
  
  
  ## write dataframe of taxa removed
  
  n_contam_taxa <- rep(n_contam_taxa, each=88*3)
  sample_cont_overlap <- rep(sample_cont_overlap, each=3)
  input_jsd <- rep(input_jsd, each=3)
  
  
  
  taxa_venn_df <- data.frame( do.call(c, removed_taxa$SCRUB ), 
                              do.call(c, removed_taxa$dec ), 
                              do.call(c, removed_taxa$`dec (lb)` ), 
                              do.call(c, removed_taxa$microDecon ))
  
  colnames( taxa_venn_df ) <- c('SCRuB', 'Decontam', 'Decontam (LB)', 'microDecon')
  
  taxa_venn_df %>% write.csv('../results/data/Fig_1/No_contamination_removed_taxa.csv')
  
  
}

run_s2_simulations(samples, controls, n_experiments=50)





noise_testing_simulations <- function(full_samps, 
                                      full_conts, 
                                      use_uploaded_samples=T,
                                      n_experiments=1){
  
  total_experiments <- n_experiments * 3  * 4
  
  n_taxa <- ncol(full_samps)
  n_samples <- 88
  n_contam_species <- NA
  
  sitch_results <- NA
  idx <- 1
  for( noise_lvl in c('low', 'med', 'high')){
    control_noise <- noise_lvl
    sample_noise <- noise_lvl
    contam_in_sample_noise <- noise_lvl
    
  for(contam_level in c(.1) ){

    for(experiment in 1:n_experiments){
      
      n_controls <- 8
      
      
      samples <- get_samps_from_upload(full_samps, n_samples) %>% get_rescaled_mat
      row.names(samples) <- paste0('Sample_', 1:n_samples)
      
      true_contamination <- full_conts[sample(nrow(full_conts), 1),] %>% rescale
      n_taxa <- ncol(samples)
      
      
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
                                 contam_in_sample_noise=contam_in_sample_noise) %>% return() }  
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
      
      
      for(plate_format in c('F')){ 
        well_positions <- Simulate_Plate_Positions(n_samples, 8, n_wells=8, plate_format=plate_format)
        well_dists <- dist(well_positions) %>% as.matrix()
        for(well_level in c(0.05)){
          well_species_limit <- NA
          observed_controls <- get_final_conts(controls_pre_well_to_well,
                                               samples,
                                               well_dists, 
                                               well_level, 
                                               well_species_limit)
          colnames(observed_samples) <- colnames(samples)
          colnames(observed_controls) <- colnames(samples)
          for(n_controls in c(8, 4, 2, 1)){ # test across 1:8 controls
            if(n_controls>1){ out <- run_SCRUB_DEC_and_RESTR(observed_samples, observed_controls[1:n_controls,], well_dists, dist_threshold=1.5)
            }else{
              tmp <- observed_controls[1, ] %>% t()
              row.names(tmp) <- c('Control_1')
              out <- run_SCRUB_DEC_and_RESTR(observed_samples, tmp, well_dists, dist_threshold=1.5)
            }
            # store_finalized_results_using_jsd
            print(paste('Running Fig_1 iteration', 
                        idx, 'of', total_experiments))
            sitch_results <- store_f1_results(sitch_results, 
                                              out,
                                              samples, 
                                              observed_samples, 
                                              n_controls, 
                                              well_species_limit, 
                                              plate_format, 
                                              well_level, 
                                              contam_level, 
                                              true_contamination,
                                              noise_lvl,
                                              idx)
            
            idx <- idx+1
          }
        }
      }
      write.csv(sitch_results, file='../results/data/Fig_1/Noise_varying_simulations.csv')
      print(idx)
    }
  }
  }
  return(sitch_results)
}





# noise_tests <-  noise_testing_simulations(samples, 
#                                           controls,
#                                           n_experiments=10)


