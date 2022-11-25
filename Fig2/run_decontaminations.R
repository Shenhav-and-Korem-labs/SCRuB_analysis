library(tidyverse)
library(magrittr)
library(glmnet)
library(torch)
install_torch()
library(decontam)
library(microDecon)
source('SCRuB/lsq_initializations.R')
source('SCRuB/spatial_functions.R')
source('SCRuB/main_functions.R')
source('Fig2/helper_functions.R')

CLEAN_SAMPLES_MICRODECON <- function(smps, cnts){
  tmp <- data.frame( rbind( colnames(smps), cnts, smps ) %>% t() )
  tmp[, 2:( ncol(tmp) ) ]  <- as.numeric( as.matrix( tmp[, 2:( ncol(tmp) ) ] ) )
  decontaminated <- decon(data = tmp, numb.blanks=nrow(cnts), numb.ind= c(nrow(smps)), taxa=F)
  
  md_out <- smps*0
  md_out[,decontaminated$decon.table$V1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t()
  
  return(list(decontaminated_samples=md_out))
}


# Reading the data
tbl <- read.csv('../data/Fig2/ASVs.csv', row.names=1) %>% as.matrix()

q <- (read.csv('../data/Fig2/SraRunTable.txt', row.names=1) %>% select(  anonymized_name ) )[row.names(tbl),] %>% as.character
row.names(tbl) <- trimws(q) %>% str_replace_all(' ', '_') 


taxa <- read.csv('../data/Fig2/ASV_alternate_taxonomy.csv', row.names = 1) %>% as.matrix
row.names(taxa) <- NULL
tax <- taxa %>% apply(MARGIN=1, function(x) paste(x, collapse=' '))
colnames(tbl) <- tax


neg <- tbl[ row.names(tbl) %>% str_detect('extractNTC'), ]
samps <- tbl[ ( row.names(tbl) %>% str_detect('extractNTC') ) == FALSE, ]

# set a seed
set.seed(1)
plate='P1'
P4 <- tbl[ which( tbl %>% row.names() %>% str_detect( plate ) ), ] 

P4 <- P4[, which(apply(P4, MARGIN=2, sum) > 0) ]
P4 <- P4[which(apply(P4, MARGIN=1, sum) > 500), ]
tax <- tax[P4 %>% apply(MARGIN=2, sum) %>% order(decreasing = TRUE)]
P4 <- P4[,P4 %>% apply(MARGIN=2, sum) %>% order(decreasing = TRUE)]

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

well_dists <- positions %>% dist(method='euclidean') %>% as.matrix()


SCR <- set_up_spatial_decontamination(P4, well_dists)

P4 <- P4[ which(apply(P4, MARGIN = 1, sum) > 500 ), ]

positions <- get_positions_df(P4)

decontam_out <- CLEAN_SAMPLES_DECONTAM(samples =  P4[row.names(P4) %>% str_detect('extractNTC')==FALSE,] %>% unname, 
                                       controls = P4[row.names(P4) %>% str_detect('extractNTC'),] %>% unname
)

decontam_lowbiomass_out <- CLEAN_SAMPLES_DECONTAM_low_biomass(samples =  P4[row.names(P4) %>% str_detect('extractNTC')==FALSE,] %>% unname, 
                                       controls = P4[row.names(P4) %>% str_detect('extractNTC'),] %>% unname
)

md_setup <- P4
colnames(md_setup) <- paste0('OTU_', 1:ncol(md_setup))
microdec_out <- CLEAN_SAMPLES_MICRODECON(smps =  md_setup[row.names(md_setup) %>% str_detect('extractNTC')==FALSE,], 
                                         cnts = md_setup[row.names(md_setup) %>% str_detect('extractNTC'),])
colnames(microdec_out$decontaminated_samples) <- colnames(P4)


REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X = samples
  to_remove <- which( (controls %>% colSums() )  > 0  )
  X[,to_remove] <- 0
  return(X)
}

restrictive_out <- REMOVE_ALL_CONTROL_SPECS( samples =  P4[row.names(P4) %>% str_detect('extractNTC')==FALSE,] %>% unname, 
                                             controls = P4[row.names(P4) %>% str_detect('extractNTC'),] %>% unname
                                             )


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



positions_ <- data.frame(wells, indices) %>%
  mutate(indices=as.integer(indices)) %>% 
  as.matrix()  

row.names(positions_) <- row.names(P4)

q_well_dists <- positions_ %>% dist(method='euclidean') %>% as.matrix()
full_set_out <- spatial_SCRUB(P4, 
                                       row.names(P4) %>% str_detect('extractNTC'), 
                                       q_well_dists, 
                                       dist_threshold = 1.5, 
                                       a_init = .95, 
                                       print_loglikelihood = F
)

row.names(full_set_out$decontaminated_samples) <- row.names(P4[row.names(P4) %>% str_detect('extractNTC') == F,])

dir.create('../results/data/Fig_2')

P4 %>% write.csv('../results/data/Fig_2/Raw.csv')
full_set_out$decontaminated_samples %>% write.csv('../results/data/Fig_2/SCRUB.csv')
decontam_out$decontaminated_samples %>% write.csv('../results/data/Fig_2/Decontam.csv')
decontam_lowbiomass_out$decontaminated_samples %>% write.csv('../results/data/Fig_2/Decontam_LB.csv')
restrictive_out %>% write.csv('../results/data/Fig_2/Restrictive.csv')
microdec_out$decontaminated_samples %>% write.csv('../results/data/Fig_2/microDecon.csv')


full_set_out$gamma %>% write.csv('../results/data/Fig_2/scrub_contaminant.csv')
1-decontam_out$cont_idx$p %>% write.csv('../results/data/Fig_2/decont_cont.csv')
decontam_lowbiomass_out$cont_idx$p %>% write.csv('../results/data/Fig_2/decontam_low_biomass_cont.csv')

(1 - microdec_out$decontaminated_samples /  md_setup[row.names(md_setup) %>% str_detect('extractNTC')==FALSE,] ) %>% colMeans(na.rm = T) %>% unname() %>%
  write.csv('../results/data/Fig_2/microdecon_cont.csv')



