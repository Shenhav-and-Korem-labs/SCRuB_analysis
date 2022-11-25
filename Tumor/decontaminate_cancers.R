

library(tidyverse)
# library(phyloseq)
library(decontam)
library(vegan)
library(ggvenn)
library(glmnet)
library(microDecon)

source('SCRuB/lsq_initializations.R')
source('SCRuB/main_functions.R')
source('SCRuB/spatial_functions.R')

set.seed(1)

load('../data/Fig3/cancer_decont_input.Rdata')

# define names of tissues
control_tissues <- c( 'control', 'NTC', 'paraf control' )
sample_tissues <- c( 'Bone', 'Breast','Colon', 'GBM', 'Lung','Melanoma', 'Ovary', 'Pancreas')


# generic functions to run facilitate running FEAST in different orders

decontam_group <- function(group, control_type, clean_func=SCRUB){
  tissues <- group %>% pull(tissue.type)
  control_tissues <- c( 'control', 'NTC', 'paraf control' )
  sample_tissues <- c( 'Bone', 'Breast','Colon', 'GBM', 'Lung','Melanoma', 'Ovary', 'Pancreas')
  #make sure the grup has the controls/tissues we need to run FEAST
  if( ( ( c( sample_tissues) %in% tissues ) %>% sum() ) * ( control_type %in% tissues ) ){ 
    getwd
    if(control_type == 'paraf control'){
      contams <- group %>% 
        filter(tissue.type == control_type) %>%
        sample_n(3) %>%
        select(-new_SEQ_NAMES, 
               -tissue.type,
               -Center,
               -DNA.Extraction.Batch,
               -PCR...New.Batch
        ) %>% as.matrix()
    }else{
      contams <- group %>% 
        filter(tissue.type ==control_type) %>%
        select(-new_SEQ_NAMES, 
               -tissue.type,
               -Center,
               -DNA.Extraction.Batch,
               -PCR...New.Batch
        ) %>% as.matrix()
    }
    
    
    
    samps <- group %>% filter(tissue.type %in% sample_tissues) %>%
      select(-new_SEQ_NAMES, 
             -tissue.type,
             -Center,
             -DNA.Extraction.Batch,
             -PCR...New.Batch
      ) %>% as.matrix()
    
    preds <- clean_func(samps, contams)
    # CLEAN_SAMPLES_lsq_init(samps, contams, new_feast = TRUE, use_adam = FALSE, COVERAGE = 10000, adam_iters = 3000)
    
    # preds <- CLEAN_SAMPLES_DECONTAM(samps, contams) #this line to run the decontam method
    return(preds)
  }
  return( 
    list(alpha = NA, 
         decontaminated_samples =  group %>%    
           filter(tissue.type !=control_type) %>%
           select(-new_SEQ_NAMES, 
                  -tissue.type,
                  -Center,
                  -DNA.Extraction.Batch,
                  -PCR...New.Batch
           ) %>% as.matrix(), 
         esimated_sources = group %>%    
           filter(tissue.type == control_type) %>%
           select(-new_SEQ_NAMES, 
                  -tissue.type,
                  -Center,
                  -DNA.Extraction.Batch,
                  -PCR...New.Batch
           ) %>% as.matrix()) )
}



recombine_into_df <- function(split_groups, out_list, control_type){
  recombined <- list()
  
  for( i in 1:length(split_groups)){
    inp <-split_groups[[i]]
    out <- out_list[[i]]
    
    # match conditions we required to run feast for the group
    if(( ((inp %>% pull(tissue.type) %in% sample_tissues) %>% sum() ) >0 )*( ((inp %>% pull(tissue.type) == control_type) %>% sum() ) > 0 )){
      inp[inp$tissue.type %in% sample_tissues, 6:ncol(inp)] <- out$decontaminated_samples
      
      if(control_type=='paraf control'){
        inp[ (inp$tissue.type == control_type )[1:3], 6:ncol(inp)] <- out$estimated_sources
      }else{
        inp[inp$tissue.type == control_type, 6:ncol(inp)] <- out$estimated_sources
      }
    }
    
    
    recombined[[i]] <- inp
  }
  
  return( bind_rows(recombined) )
}




run_decontam_layer <- function(df,
                               decontam_func, 
                               control_tissue, 
                               group_level, 
                               parallelize = TRUE){
  group_level <- as.symbol(group_level)
  split_groups <-  df %>% split( df %>% pull(group_level))
  
  if(parallelize == TRUE){
    numCores <- detectCores()
    cl <- makeCluster(numCores)
    
    clusterEvalQ(cl, {
      library(tidyverse)
      source('../beta_functions.R')
      source('../Well_to_well/spatial_functions.R')
    })
    
    outputs <- parLapply(cl,
                               X = split_groups, 
                               fun = decontam_group, 
                               control_type = control_tissue
    )
  }else{ outputs <- lapply( split_groups, decontam_group, control_type = control_tissue, clean_func=decontam_func) }
  
  recombined_df <- recombine_into_df(split_groups, outputs, control_tissue)
  return(recombined_df)
}

set.seed(2)
reverse_paraf_out <- all_info %>%
  run_decontam_layer( SCRUB, 'paraf control', 'Center', parallelize = FALSE)

print('did SCRUB paraf')

reverse_ntc_out <- reverse_paraf_out %>%  #all_info %>% 
  run_decontam_layer( SCRUB,  'NTC', "PCR...New.Batch", parallelize = FALSE)

print('did SCRUB NTC')

reverse_cont_out <-  reverse_ntc_out %>%
  run_decontam_layer( SCRUB,  'control', "DNA.Extraction.Batch", parallelize = FALSE)

dir.create('../results/data/Tumor')

print('finished SCRUB')

write.csv( reverse_cont_out, file = '../results/data/Tumor/SCRUB_result.csv')


their_output <- (read.csv('../data/Fig3/df_freqs.csv') ) %>% t()
their_output[is.na(their_output)] <- 0

tissue_names <- all_info %>% filter(tissue.type %in% sample_tissues) %>% pull(new_SEQ_NAMES)

row.names( their_output ) <- row.names( their_output ) %>% str_replace('X2', '2')
their_output[their_output %>% row.names() %in% pull( all_info, new_SEQ_NAMES), ] %>% nrow()
filter_tracker <- read.csv('../data/Fig3/df_filter_tracker.csv')

their_output <-  their_output[, filter_tracker %>% pull(species)!= '' ]

all_inf_phylos <- read.csv('../data/Fig3/df_freqs_before_global_contaminants_filter.csv')[2:8]

all_inf_phylos$all_inf_idx <- seq(1, nrow(all_inf_phylos))
filter_tracker_phylos <- filter(filter_tracker, species != '')[, 2:8]
filter_tracker_phylos$filt_idx <- seq(1, nrow(filter_tracker_phylos))

merged_idx <- all_inf_phylos %>%
  merge(filter_tracker_phylos, 
        by.x= c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
        by.y= c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
        all=TRUE) %>% 
  select( all_inf_idx, filt_idx ) %>%
  arrange(all_inf_idx)

merged_idx$filt_idx[ merged_idx$filt_idx %>% is.na() ] <- seq( 1 + nrow(filter_tracker_phylos), nrow(all_inf_phylos) )


num_pads <- nrow(all_inf_phylos) - nrow(filter_tracker_phylos)
their_out_padded <- their_output %>% cbind( matrix( 0,
                                                    their_output %>% nrow(), 
                                                    num_pads) )

matching_format <- their_out_padded[2:nrow(their_out_padded), merged_idx$filt_idx]


#use this cell to run decontam
decontam_paraf_out <- all_info %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM,  'paraf control', 'Center', parallelize = FALSE)

print('did decontam paraf')

decontam_ntc_out <- decontam_paraf_out %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM,  'NTC', "PCR...New.Batch", parallelize = FALSE)

print('did decontam NTC')

decontam_cont_out <-  decontam_ntc_out %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM, 'control', "DNA.Extraction.Batch", parallelize = FALSE)

print('finished decontam')


write.csv(decontam_cont_out, file = '../results/data/Tumor/decontam_result.csv')

REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X <- samples
  to_remove <- which( (controls %>% colSums() )  >  0  )
  X[,to_remove] <- 0
  
  return(list(decontaminated_samples=X, 
              estimated_sources=colSums(controls)))
}

neg_idx <- all_info$tissue.type %in% c('control', 'NTC', 'paraf control') 
global_fully_restrictive <- cbind(all_info[neg_idx==F, 1:5],
                                  REMOVE_ALL_CONTROL_SPECS(all_info[neg_idx==F, 6:ncol(all_info)], 
                                                           all_info[neg_idx, 6:ncol(all_info)])$decontaminated_samples
                                  )
  
global_fully_restrictive %>% write.csv('../results/data/Tumor/global_fully_restrictive_result.csv')

CLEAN_SAMPLES_MICRODECON <- function(smps, cnts){
  cname_placeholder <- colnames(smps)
  colnames(smps) <- paste0('OTU_', 1:ncol(smps))
  colnames(cnts) <- paste0('OTU_', 1:ncol(smps))
  print('new run')
  print(sum(smps))
  
  if(sum(smps)==0){return(list(decontaminated_samples=smps, 
                               estimated_sources=colSums(cnts))) }
  tmp <- data.frame( rbind( colnames(smps), cnts, smps ) %>% t() )
  tmp[, 2:( ncol(tmp) ) ]  <- as.numeric( as.matrix( tmp[, 2:( ncol(tmp) ) ] ) )
  decontaminated <- decon(data = tmp, numb.blanks=nrow(cnts), numb.ind= c(nrow(smps)), taxa=F)
  # print( as.matrix(decontaminated$decon.table)[1:20, ] )
  md_out <- smps*0
  md_out[,decontaminated$decon.table$X1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t() %>% as.double()
  print(dim(decontaminated$decon.table))
  print(sum(md_out))
  colnames(md_out) <- cname_placeholder
  return(list(decontaminated_samples=md_out, 
              estimated_sources=colSums(cnts)))
}

set.seed(3)
microdecon_paraf_out <- all_info %>%
  run_decontam_layer( CLEAN_SAMPLES_MICRODECON,  'paraf control', 'Center', parallelize = FALSE)

microdecon_ntc_out <- microdecon_paraf_out %>%
  run_decontam_layer( CLEAN_SAMPLES_MICRODECON,  'NTC', "PCR...New.Batch", parallelize = FALSE)

microdecon_control_out <- microdecon_ntc_out %>%
  run_decontam_layer( CLEAN_SAMPLES_MICRODECON,  'control', "DNA.Extraction.Batch", parallelize = FALSE)


microdecon_control_out %>% write.csv(file = '../results/data/Tumor/microdecon_result.csv')


set.seed(100)
#use this cell to run decontam
decontam_paraf_out_lb <- all_info %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM_low_biomass,  'paraf control', 'Center', parallelize = FALSE)

print('did decontam paraf')

decontam_ntc_out_lb <- decontam_paraf_out_lb %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM_low_biomass,  'NTC', "PCR...New.Batch", parallelize = FALSE)

print('did decontam NTC')

decontam_cont_out_lb <-  decontam_ntc_out_lb %>%
  run_decontam_layer( CLEAN_SAMPLES_DECONTAM_low_biomass,  'control', "DNA.Extraction.Batch", parallelize = FALSE)


write.csv(decontam_cont_out_lb, file = '../results/data/Tumor/decontam_low_biomass_result.csv')


melanoma_idcs <- which( ( ( all_info %>%
                              filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)]  ) %>% as.matrix() %>%
                          colSums() > 1e-2 ) %>% unname()


div_ind <- 'shannon'

diver_raw <- ( ( all_info %>%
                   filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)]  ) %>% as.matrix() %>%
  diversity(index=div_ind)

diver_scrub <- ( ( reverse_cont_out %>%
                     filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)]  ) %>% as.matrix() %>%
  # specnumber()
  diversity(index=div_ind)

diver_rest <- matching_format[which( row.names(matching_format) %in% 
                                       ( all_info %>% filter( tissue.type == 'Melanoma') )$new_SEQ_NAMES ), ] %>%
  diversity(index=div_ind) %>% unname()

diver_dec <- ( ( decontam_cont_out %>%
                   filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
  # specnumber()
  diversity(index=div_ind)


full_restrictive_div <- ( ( global_fully_restrictive %>%
                      filter(tissue.type=='Melanoma') )[, 6:ncol(global_fully_restrictive)] ) %>%
  diversity(index=div_ind)


diver_dec_lb <- ( ( decontam_cont_out_lb %>%
                   filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
  # specnumber()
  diversity(index=div_ind)


diver_microdec <- ( ( microdecon_control_out %>%
        filter(tissue.type=='Melanoma') )[, 6:ncol(microdecon_control_out)] ) %>%
         diversity(index=div_ind)


n_melanoma_samples <- length(diver_raw)
divers_df <- data.frame( 
  c( diver_raw, diver_scrub, diver_dec,  diver_rest, diver_dec_lb, 
      full_restrictive_div, diver_microdec ), 
  c( rep('Raw', n_melanoma_samples), rep('SCRUB', n_melanoma_samples), 
     rep('Decontam', n_melanoma_samples), rep( 'Restrictive (Original)', n_melanoma_samples ),
     rep('Decontam (LB)', n_melanoma_samples), 
     rep('Restrictive', n_melanoma_samples), 
     rep('microDecon', n_melanoma_samples) )
     )

colnames(divers_df) <- c('Shannon', 'Dataset')

divers_df %>% write.csv('../results/data/Tumor/Shannon_diversities.csv')

dir.create('../results/Plots/')
dir.create('../results/Plots/fig_3/')

ggvenn(
  list(
    "Custom"= ((matching_format[which( row.names(matching_format) %in% 
                                                         ( all_info %>% filter( tissue.type == 'Melanoma') )$new_SEQ_NAMES ), ] %>%
                                  colSums )[melanoma_idcs] == 0 ) %>% which , 
    "SCRUB" = ( ( ( ( reverse_cont_out %>%
                        filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                    colSums )[melanoma_idcs] == 0 ) %>% which ,
    "Decontam (LB)"=( ( ( ( decontam_cont_out_lb %>%
                              filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                          colSums )[melanoma_idcs] == 0 ) %>% which, 
    "Restrictive"=( ( ( ( global_fully_restrictive %>%
                         filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                     colSums )[melanoma_idcs] == 0 ) %>% which
  ), 
  columns = c("Custom", 
              "SCRUB" , 
              "Decontam (LB)", 
              "Restrictive"),
   fill_color=c('#FF69B4', '#c44e52', 'darkgreen', '#dd8452'),
  show_percentage = F, 
  fill_alpha = .7
)
ggsave('../results/Plots/fig_3/Melanoma_4_group_VD_with_numbers.pdf', device='pdf', dpi=900)


ggvenn(
  list( "Custom"= 1:10, "SCRUB" =1:10,"Decontam (LB)"=1:10,"Restrictive"=1:10),
  columns = c("Restrictive (Original)" ,
                     "SCRUB" ,
                     "Decontam (LB)",
                      "Restrictive"),
  fill_color=c('#FF69B4', '#c44e52', 'darkgreen', '#dd8452'),
  show_percentage = F, 
  fill_alpha = .8, 
  show_elements=F, 
  text_size=0, 
  set_name_size=0 ) 

ggsave('../results/Plots/fig_3/Melanoma_4_group_VD.pdf', device='pdf', dpi=900)


