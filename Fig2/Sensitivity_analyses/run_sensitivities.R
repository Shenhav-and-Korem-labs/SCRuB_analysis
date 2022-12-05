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

full_sensitivity_results <- c()

all_raw_samples <- list()
all_scrubbed_samples <- list()
i <- 0
for(ind in row.names(P4[ (row.names(P4) %>% str_detect('extractNTC') )==F,])){
  i <- i+1
print(ind)
P4_copy <- P4
P4_copy[ind, ] <- P4_copy["P1.HT1PCRA.SourceEcoli.G10", ]
  
full_set_out <- spatial_SCRUB(P4_copy, 
                                       row.names(P4) %>% str_detect('extractNTC'), 
                                       q_well_dists, 
                                       dist_threshold = 1.5, 
                                       a_init = .95, 
                                       print_loglikelihood = F )


all_raw_samples[[i]] <- P4_copy[ (row.names(P4) %>% str_detect('extractNTC') )==F,]
all_scrubbed_samples[[i]] <- full_set_out$decontaminated_samples


}

# full_sensitivity_results %>% write.csv('../results/data/Fig_2/Sensitivity/Iterative_results.csv')


full_ecoli_scrubbed <- data.table::rbindlist(all_scrubbed_samples %>% lapply(function(x){
                                                                          tmp <- data.frame(x) 
                                                                          if(nrow(tmp)==0)return(tmp)
                                                                          row.names(tmp) <- rns
                                                                          return(tmp) }
                                                                        ) )
full_ecoli_raw <- data.table::rbindlist(all_raw_samples %>% lapply(function(x){
                                                                        tmp <- data.frame(x) 
                                                                        if(nrow(tmp)==0)return(tmp)
                                                                        row.names(tmp) <- rns
                                                                        return(tmp) }
                                                                      ) )

row.names(full_ecoli_raw) <-  rep(rns, 64) %>% make.names(unique=T)
row.names(full_ecoli_scrubbed) <-  rep(rns, 64) %>% make.names(unique=T)

full_ecoli_scrubbed %>% write.csv('../results/data/Fig_2/Sensitivity/Ecoli_scrub.csv')
full_ecoli_raw %>% write.csv('../results/data/Fig_2/Sensitivity/Ecoli_raw.csv')



ind_scrubbed <- c()
ind_raws <- c()

for(i in 1:100){
  print(i)
  inds <- row.names(P4[ (row.names(P4) %>% str_detect('Sourc') )==T,]) %>% sample(2)
  P4_copy <- P4
  P4_copy[ inds[1], ] <- P4[ inds[2], ]
  P4_copy[ inds[2], ] <- P4[ inds[1], ]
  
  full_set_out <- spatial_SCRUB(P4_copy, 
                                row.names(P4_copy) %>% str_detect('extractNTC'), 
                                q_well_dists, 
                                dist_threshold = 1.5, 
                                a_init = .95, 
                                print_loglikelihood = F )
  
  colnames(full_set_out$decontaminated_samples) <- colnames(P4_copy)
  rownames(full_set_out$decontaminated_samples) <-  row.names( P4_copy[ (row.names(P4_copy) %>% str_detect('extractNTC') ) == F, ])
  
  spec_val <- which( P4_copy[inds[1],] == max(P4_copy[inds[1],] ) )
  
  ind_scrubbed <- c(ind_scrubbed, get_rescaled_mat( full_set_out$decontaminated_samples )[inds[1], spec_val])
  ind_raws <- c(ind_raws, get_rescaled_mat( P4_copy )[inds[1], spec_val])
}
  

data.frame( 'composition'=c(ind_scrubbed, ind_raws), 
            'name'=c(rep('SCRuB', length(ind_scrubbed)), rep('No decontamination', length(ind_raws)))) %>%
  write.csv('../results/data/Fig_2/Sensitivity/V3_variation.csv')

















