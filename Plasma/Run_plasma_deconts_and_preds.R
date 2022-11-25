

install.packages('splitstackshape', repos='http://cran.us.r-project.org')
library(splitstackshape)
library(tidyverse)
library(biomformat) 
library(vegan)
library(glmnet)
library(torch)
library(microDecon)
# install_torch()

source('SCRuB/lsq_initializations.R')
source('SCRuB/spatial_functions.R')
source('SCRuB/main_functions.R')



for(seed in 0:9){
  set.seed(seed)
print(paste('running iteration', seed))
dec_ind <- read.csv('../data/Fig3_plasma/Metadata-Plasma-For-Decontam-With-Negative-And-Positive-Controls.csv')

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig3_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)
snmDataKrakenCFDecontamDPQC <-read.csv('../data/Fig3_plasma/Kraken-Plasma-Voom-SNM-Age-And-Sex-Data.csv', row.names=1)

data <- read_biom(biom_file = '../data/Fig3_plasma/136205_47212_analysis_Metagenomic_Woltkav011Databaseoptgenomeqiitadbsqpwoltkarep200Rep200BIOMnonebiom.biom')

taxa_names <- c()

for(i in 1:length(data$rows)) taxa_names <- c(taxa_names, data$rows[[i]]$id)

samp_names <- c()

for(i in 1:length(data$columns)) samp_names <- c(samp_names, data$columns[[i]]$id)

full_df <- matrix( unlist(data$data), byrow=TRUE, nrow=length(data$data) ) 
colnames(full_df) <- samp_names
row.names(full_df) <- taxa_names
full_df <- full_df %>% t()


metadata <- read.csv('../data/Fig3_plasma/47212_47212_analysis_mapping.txt', sep='\t')

remov_left <- function(x, n){
  substr(x, n, nchar(x))
}

unique_metadata <- metadata %>%
  filter(as.character(X.SampleID) %in% row.names(full_df)) %>%
  mutate(general_id = remov_left(as.character(X.SampleID), 7)) %>%
  group_by(general_id) %>% sample_n(1)


unique_samps <- full_df[as.character(unique_metadata$X.SampleID), ]



final_df <- full_df[which(row.names(full_df) %>% str_detect('Control') == F ), ] %>%
  rbind( full_df[metadata %>% filter(str_detect(X.SampleID, 'Control'), 
                                     #X.SampleID %>% str_detect('Control'),
                                     X.SampleID %>% str_detect('12691') ) %>% 
                   pull(X.SampleID) %>% unique() %>% as.character, ] ) 


removed_copies <- metadata %>% filter( ( str_detect(X.SampleID, 'Control') == F ) ) %>% #&(str_detect(description, 'HIV') == F ) ) %>%
  rbind( metadata %>% filter(str_detect(X.SampleID, 'Control'), 
                             X.SampleID %>% str_detect('12691') ) )

unique_metadata <- as.data.frame(unique_metadata)
row.names(unique_metadata) <- unique_metadata[, 'X.SampleID'] %>% as.character()
well_dists <- unique_metadata %>%
  mutate(well_loc= sample_well %>% substr(1,1) %>% sapply( function(x) which(LETTERS==x)[1]), 
         indices_loc = sample_well%>% substr(2,3) %>% as.integer ) %>%
  select(well_loc, indices_loc) %>%
  dist(method = 'euclidean') %>% as.matrix()


## run through SCRuB

tmp_fwd <- as.matrix( unique_samps[ unique_metadata %>% 
                                      filter(sample_type %>% 
                                               str_detect('control blank library prep') == F) %>% 
                                      pull(X.SampleID) %>% 
                                      as.character(), ] )


plasma_scrubbed <- spatial_SCRUB(data=tmp_fwd, 
                                 is_control = unique_metadata[row.names(tmp_fwd),]$sample_type=='control blank DNA extraction',
                                          well_dists = well_dists, 
                                          dist_threshold =1.5
)

row.names( plasma_scrubbed$decontaminated_samples ) <- unique_metadata[row.names(tmp_fwd),
][ unique_metadata[row.names(tmp_fwd),
]$sample_type!='control blank DNA extraction', ] %>% row.names()



next_lvl_mat <- plasma_scrubbed$decontaminated_samples %>%
  rbind(unique_samps[row.names(unique_metadata %>%
                                 filter(sample_type=='control blank library prep')), ])



plasma_scrubbed_2nd <- spatial_SCRUB(data=next_lvl_mat, 
  is_control = c( rep(F, nrow(plasma_scrubbed$decontaminated_samples) ),
                  rep(T, sum(unique_metadata$sample_type=='control blank library prep'))),
                                              well_dists = well_dists,
                                              dist_threshold =  1.5 )


row.names(plasma_scrubbed_2nd$decontaminated_samples) <- unique_metadata[row.names(plasma_scrubbed$decontaminated_samples),
                                                ][ unique_metadata[row.names(plasma_scrubbed$decontaminated_samples ),
                                                ]$sample_type!='control blank library prep', ] %>% row.names()


colnames(plasma_scrubbed_2nd$decontaminated_samples) <- colnames(unique_samps)

scrub_df <- plasma_scrubbed_2nd$decontaminated_samples[,  ( ( plasma_scrubbed_2nd$decontaminated_samples %>% colSums() ) > 0 ) %>% which]
scrub_df <- scrub_df[row.names(metadataPSMatchedDPQCFiltered), ]


scrub_df <- scrub_df[, (colSums(scrub_df) > 500) %>% which]


## run through microdecon

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
  md_out[,decontaminated$decon.table$V1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t() %>% as.double()
  print(dim(decontaminated$decon.table))
  print(sum(md_out))
  colnames(md_out) <- cname_placeholder
  return(list(decontaminated_samples=md_out, 
              estimated_sources=colSums(cnts)))
}

tmp_fwd <- as.matrix( unique_samps[ unique_metadata %>% 
                                      filter(sample_type %>% 
                                               str_detect('control blank library prep') == F) %>% 
                                      pull(X.SampleID) %>% 
                                      as.character(), ] )


is_control <-  unique_metadata[row.names(tmp_fwd),]$sample_type=='control blank DNA extraction'
plasma_microdec <- CLEAN_SAMPLES_MICRODECON(tmp_fwd[is_control==F, ], tmp_fwd[is_control, ] )

row.names( plasma_microdec$decontaminated_samples ) <- unique_metadata[row.names(tmp_fwd),
][ unique_metadata[row.names(tmp_fwd),
]$sample_type!='control blank DNA extraction', ] %>% row.names()



next_lvl_mat <- plasma_microdec$decontaminated_samples %>%
  rbind(unique_samps[row.names(unique_metadata %>%
                                 filter(sample_type=='control blank library prep')), ])


is_control = c( rep(F, nrow(plasma_microdec$decontaminated_samples) ),
                rep(T, sum(unique_metadata$sample_type=='control blank library prep')))

plasma_microdec_2nd <- CLEAN_SAMPLES_MICRODECON(next_lvl_mat[is_control==F, ], next_lvl_mat[is_control, ] )


row.names(plasma_microdec_2nd$decontaminated_samples) <- unique_metadata[row.names(plasma_microdec$decontaminated_samples),
][ unique_metadata[row.names(plasma_microdec$decontaminated_samples ),
]$sample_type!='control blank library prep', ] %>% row.names()


colnames(plasma_microdec_2nd$decontaminated_samples) <- colnames(unique_samps)

microdec_df <- plasma_microdec_2nd$decontaminated_samples[,  ( ( plasma_microdec_2nd$decontaminated_samples %>% colSums() ) > 0 ) %>% which]
microdec_df <- microdec_df[row.names(metadataPSMatchedDPQCFiltered), ]


microdec_df <- microdec_df[, (colSums(microdec_df) > 500) %>% which]


# Decontam, Restrictive
library(decontam)

dec_setup <- unique_samps[dec_ind$X,]
dim(dec_setup)

restrictive <- dec_setup[ , which( dec_setup[ dec_ind$decontam_prevalence_use, ] %>% colSums() == 0 ) ]
restrictive <- restrictive[, which( restrictive %>% colSums() > 500 ) ]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



dec_out <- isContaminant(seqtab = dec_setup, # final_df, 
                         neg = dec_ind$decontam_prevalence_use,#final_df %>% row.names() %>% str_detect('Control'), # final_df%>% row.names() %>% str_detect('Control'),
                         conc = removed_copies$well_conc + 1e-16,
                         method='combined',
                         threshold = .5, 
                         detailed=T, 
)

decontammed_data <- full_df[, (dec_out$contaminant==F) %>% which]
decontammed_data <- decontammed_data[row.names(metadataPSMatchedDPQCFiltered), ]
decontammed_data <- decontammed_data[ , (decontammed_data %>% colSums() > 500 ) %>% which ]


dec_out <- isContaminant(seqtab = dec_setup, # final_df, 
                         neg = dec_ind$decontam_prevalence_use,#final_df %>% row.names() %>% str_detect('Control'), # final_df%>% row.names() %>% str_detect('Control'),
                         conc = removed_copies$well_conc,
                         method='prevalence',
                         threshold = .1, 
                         detailed=T, 
)

decontammed_data_standard <- full_df[, (dec_out$contaminant==F) %>% which]
decontammed_data_standard <- decontammed_data_standard[row.names(metadataPSMatchedDPQCFiltered), ]
decontammed_data_standard <- decontammed_data_standard[ , (decontammed_data_standard %>% colSums() > 500 ) %>% which ]




dec_out <- isNotContaminant(seqtab = dec_setup, # final_df, 
                            neg = dec_ind$decontam_prevalence_use,#final_df %>% row.names() %>% str_detect('Control'), # final_df%>% row.names() %>% 
                            threshold = .5, 
                            detailed=T, 
)

decontammed_data_low_bm <- full_df[, (dec_out$not.contaminant==T) %>% which]
decontammed_data_low_bm <- decontammed_data_low_bm[row.names(metadataPSMatchedDPQCFiltered), ]
decontammed_data_low_bm <- decontammed_data_low_bm[ , (decontammed_data_low_bm %>% colSums() > 500 ) %>% which ]




vsnm <- function(qcData){
  ## Load packages ##
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  require(doMC)
  require(tibble)
  require(gbm)
  
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)
  
  qcMetadata <- metadataPSMatchedDPQCFiltered # ADAPT THIS AS NEEDED
  # qcData <- qqq # ADAPT THIS AS NEEDED
  
  # Set up design matrix
  covDesignNorm <- model.matrix(~0 + disease_type_consol + #randomized_disease_type_consol + #disease_type_consol +
                                  host_age + # host_age should be numeric
                                  sex, # sex should be a factor
                                data = qcMetadata)
  
  # Check row dimensions
  dim(covDesignNorm)[1] == dim(qcData)[1]
  
  # The following corrects for column names that are incompatible with downstream processing
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Quantile normalize and plug into voom
  dge <- DGEList(counts = counts)
  vdge <<- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, 
                normalize.method="quantile")
  
  # List biological and normalization variables in model matrices
  bio.var <- model.matrix(~disease_type_consol, #randomized_disease_type_consol, #disease_type_consol,
                          data=qcMetadata)
  
  adj.var <- model.matrix(~host_age +
                            sex,
                          data=qcMetadata)
  
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  print(dim(adj.var))
  print(dim(bio.var))
  print(dim(t(vdge$E)))
  print(dim(covDesignNorm))
  
  snmDataObjOnly <- snm(raw.dat = vdge$E, 
                        bio.var = bio.var, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE,
                        diagnose = TRUE)
  snmData <<- t(snmDataObjOnly$norm.dat)
  
}

scrubbed_normalized <- vsnm(scrub_df[row.names(metadataPSMatchedDPQCFiltered), ])

microdecon_normalized <- vsnm(microdec_df[row.names(metadataPSMatchedDPQCFiltered), ])

# scrubbed_normalized <- vsnm(scrub_df[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ), ] )
# 
raw_inp <- unique_samps[row.names(metadataPSMatchedDPQCFiltered),]
raw_inp <- raw_inp[, colSums(raw_inp)>500]
# raw_inp <- group_to_genus(raw_inp)
raw_normalized <- vsnm(raw_inp)

dec_normalized <- vsnm(decontammed_data[row.names(metadataPSMatchedDPQCFiltered), ])#
# dec_normalized <- vsnm(decontammed_data[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])


dec_standard_normalized <- vsnm(decontammed_data_standard[row.names(metadataPSMatchedDPQCFiltered), ])#

# dec_standard_normalized <- vsnm(decontammed_data_standard[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])

dec_lb_normalized <- vsnm(decontammed_data_low_bm[row.names(metadataPSMatchedDPQCFiltered), ])#
# dec_lb_normalized <- vsnm(decontammed_data_low_bm[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])



restrictive_normalized <-  vsnm( restrictive[row.names(metadataPSMatchedDPQCFiltered) %>% as.character(), ])#


## CODE FROM KNIGHT LAB'S GITUHB


# Load dependencies
require(devtools)
require(doMC)
require(tibble)
require(gbm)
require(splitstackshape)
require(reshape2)
require(ggpubr)
require(caret) # for model building
require(pROC) # for AUC calculations
require(purrr) # for functional programming using map()
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(cowplot) # for plotting
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(caret) # for machine learning

defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)
customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                              n.trees = floor((1:3) * 50),
                              shrinkage = 0.1,
                              n.minobsinnode = 1)
numKFold <- 4
numResampleIter <- 1

ml2DTs <- function(snmData, 
                   classOfInterest = "Lung Adenocarcinoma", 
                   cutPoint = 0.5, 
                   samplingSize = 20, 
                   caretTuneGrid = defaultGBMGrid){
  
  metaTmp1 <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% c("PRAD",
                                                                                                                 "SKCM",
                                                                                                                 "NSCLC")),])
  tmp <- metaTmp1
  tmp$disease_type_consol <- factor(ifelse(metaTmp1$disease_type_consol == classOfInterest, yes = classOfInterest, no = "Other"))
  metadataSimSampled <- as.data.frame(stratified(tmp,
                                                 group = "disease_type_consol",
                                                 size = samplingSize,
                                                 keep.rownames = TRUE,
                                                 replace = FALSE,
                                                 bothSets = FALSE))
  rownames(metadataSimSampled) <- metadataSimSampled$rn
  mlDataY <- metadataSimSampled
  mlDataX <- snmData[rownames(mlDataY),]
  
  set.seed(seed)
  index <- createDataPartition(mlDataY$disease_type_consol, p = 0.7, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$disease_type_consol
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$disease_type_consol
  # print(testY)
  
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
  
  set.seed(seed)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = "gbm",
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = customGBMGrid)
  
  positiveClass <- gsub(" ","", classOfInterest)
  negativeClass <- "Other"
  
  predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
  fg <- predProbs[refactoredTestY == positiveClass]
  bg <- predProbs[refactoredTestY == negativeClass]
  
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  # par(mfrow = c(1,2))
  plot(prroc_roc)
  plot(prroc_pr)
  # dev.off()
  
  
  predClass <- predict(mlModel, newdata = testX)
  
  confusionMatrix(table(predict(mlModel, newdata = testX, type="prob")[,positiveClass] >= cutPoint,
                        refactoredTestY == positiveClass))
}

#-----------------------------------------#
# Machine learning
#-----------------------------------------#


# mlHvsC <- function(snmData){
# Load dependencies


mlHvsC <- function(snmData){
  
  numCores <- detectCores()
  registerDoMC(cores=numCores)
  
  defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                 n.trees = floor((1:3) * 50),
                                 shrinkage = 0.1,
                                 n.minobsinnode = 5)
  customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                n.trees = floor((1:3) * 50),
                                shrinkage = 0.1,
                                n.minobsinnode = 1)
  
  caretTuneGrid <- defaultGBMGrid
  numKFold <- 4
  numResampleIter <- 1
  
  mlDataY <- metadataPSMatchedDPQCFiltered
  mlDataX <- snmData[rownames(mlDataY),]
  
  set.seed(seed)
  index <- createDataPartition(mlDataY$HvsC, p = 0.7, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$HvsC
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$HvsC
  
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
  
  set.seed(seed)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = "gbm",
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = defaultGBMGrid)
  
  positiveClass <- "Cancer"
  negativeClass <- "Control"
  predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
  fg <- predProbs[refactoredTestY == positiveClass]
  bg <- predProbs[refactoredTestY == negativeClass]
  
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  plot(prroc_roc)
  plot(prroc_pr)
  
  predClass <- predict(mlModel, newdata = testX)
  print(confusionMatrix(data = predClass, reference = refactoredTestY, positive = positiveClass))
}

#-----------------------------------------------------

loocvDTs <- function(snmData, samplingSize = 15, DTs, caretTuneGrid = defaultGBMGrid,
                     filenameString = paste(DTs,collapse = "__"), HvsCFlag = FALSE){
  
  if(HvsCFlag){
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    metaTmpX$disease_type_consol <- metaTmpX$HvsC
    classes <- gsub(" ","",levels(metaTmpX$disease_type_consol))
  } else{
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    classes <- gsub(" ","",DTs)
  }
  
  # Do LOOCV model building and testing
  
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  metaData <- metaTmpX
  snmData <- snmData # dataPSUniqueDecontamQC # 
  iterSize <- 1
  for(jj in 1:iterSize){
    metadataSimSampled <- as.data.frame(stratified(metaData,
                                                   group = "disease_type_consol",
                                                   size = samplingSize,
                                                   keep.rownames = TRUE,
                                                   replace = FALSE,
                                                   bothSets = FALSE))
    rownames(metadataSimSampled) <- metadataSimSampled$rn
    mlDataY <- metadataSimSampled
    mlDataX <- snmData[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    # Create data partitions
    # set.seed(42)
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    obsClass <- vector()
    predClass <- vector()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    for(ii in 1:length(indexSuper)){
      print(sprintf("Iteration: %d/%d", ii, length(indexSuper)))
      index <- indexSuper[ii]
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$disease_type_consol
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$disease_type_consol
      
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
      
      obsClass[ii] <- as.character(refactoredTestY)
      
      set.seed(seed)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = multiClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = FALSE,
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)
      
      predProbs[ii] <- list(predict(mlModel, newdata = testX, type = "prob"))
      predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))
      
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      varImpBestModelDF2OrderedNonzeroList[[ii]] <- varImpBestModelDF2OrderedNonzero
      
      rm(mlModel)
    }
    
    loocvPreds <- cbind(obs = factor(obsClass,
                                     levels = classes),
                        pred = factor(predClass,
                                      levels = classes),
                        do.call(rbind,predProbs))

    
    multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)
    print(multiClassSummaryStats[[jj]])
    
    loocvPreds %>% write.csv( paste0(filenameString, "__Preds.csv"))
    
    filenameROC <- paste0(filenameString,"__ROC.png")
    filenamePR <- paste0(filenameString,"__PR.png")
    filenameROCData <- paste0(filenameString,"__Data__ROC.csv")
    filenamePRData <- paste0(filenameString,"__Data__PR.csv")
    filenameSink <- paste0(filenameString,"__CM.txt")
    
    predProbs <- loocvPreds[,DTs[1]]
    fg <- predProbs[loocvPreds$obs == DTs[1]]
    bg <- predProbs[loocvPreds$obs == DTs[2]]
    
    prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    
    png(filename=filenameROC, width = 6, height = 4, units = 'in', res = 300)
    plot(prroc_roc)
    dev.off()
    
    png(filename=filenamePR, width = 6, height = 4, units = 'in', res = 300)
    plot(prroc_pr)
    dev.off()
    
    rocCurveData <- cbind(as.data.frame(prroc_roc$curve), DT1 = DTs[1], DT2 = DTs[2])
    prCurveData <- cbind(as.data.frame(prroc_pr$curve), DT1 = DTs[1], DT2 = DTs[2])
    
    write.table(prCurveData, sep=",", file = filenamePRData, col.names = FALSE)
    write.table(rocCurveData, sep=",", file = filenameROCData, col.names = FALSE)
  }
  
  print(confusionMatrix(loocvPreds$obs, loocvPreds$pred))
  multiClassSummaryStatsDist <- data.frame(do.call(rbind, multiClassSummaryStats))
  
  sink(filenameSink)
  print(print(confusionMatrix(loocvPreds$obs, loocvPreds$pred)))
  sink()
  
  return(multiClassSummaryStats[jj])
}

#-----------------------------------------------------

ml2DTs <- function(snmData, 
                   classOfInterest = "Lung Adenocarcinoma", 
                   cutPoint = 0.5, 
                   samplingSize = 20, 
                   caretTuneGrid = defaultGBMGrid){
  
  metaTmp1 <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% c("PRAD",
                                                                                                                 "SKCM",
                                                                                                                 "NSCLC")),])
  tmp <- metaTmp1
  tmp$disease_type_consol <- factor(ifelse(metaTmp1$disease_type_consol == classOfInterest, yes = classOfInterest, no = "Other"))
  metadataSimSampled <- as.data.frame(stratified(tmp,
                                                 group = "disease_type_consol",
                                                 size = samplingSize,
                                                 keep.rownames = TRUE,
                                                 replace = FALSE,
                                                 bothSets = FALSE))
  rownames(metadataSimSampled) <- metadataSimSampled$rn
  mlDataY <- metadataSimSampled
  mlDataX <- snmData[rownames(mlDataY),]
  
  set.seed(seed)
  index <- createDataPartition(mlDataY$disease_type_consol, p = 0.7, list = FALSE)
  trainX <- mlDataX[index,]
  trainY <- mlDataY[index,]$disease_type_consol
  testX <- mlDataX[-index,]
  testY <- mlDataY[-index,]$disease_type_consol

  
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
  
  set.seed(seed)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = "gbm",
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   verbose = TRUE,
                   metric = "ROC",
                   tuneGrid = customGBMGrid)
  
  positiveClass <- gsub(" ","", classOfInterest)
  negativeClass <- "Other"
  
  predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
  fg <- predProbs[refactoredTestY == positiveClass]
  bg <- predProbs[refactoredTestY == negativeClass]
  
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  # par(mfrow = c(1,2))
  plot(prroc_roc)
  plot(prroc_pr)
  # dev.off()
  
  
  predClass <- predict(mlModel, newdata = testX)
  
  confusionMatrix(table(predict(mlModel, newdata = testX, type="prob")[,positiveClass] >= cutPoint,
                        refactoredTestY == positiveClass))
}

dir.create('../results/data/Fig3_plasma')

# Shannon diversities

div_ind <- 'shannon'

n_plasma_samples <- nrow(snmDataKrakenCFDecontamDPQC)

diver_raw <-raw_inp %>%
  diversity(index=div_ind)

diver_scrub <- rbind(scrub_df) %>% cbind( matrix(0,nrow=nrow(scrub_df), ncol = ncol(raw_inp)-ncol(scrub_df) ) ) %>%
  diversity(index=div_ind)

diver_rest <- restrictive[row.names(raw_inp),]  %>%
  diversity(index=div_ind)


diver_dec <- rbind(decontammed_data) %>% cbind( matrix(0,nrow=nrow(decontammed_data), ncol = ncol(raw_inp)-ncol(decontammed_data) ) ) %>%
  diversity(index=div_ind)


diver_lb_dec <- rbind(decontammed_data_low_bm) %>% cbind( matrix(0,nrow=nrow(decontammed_data_low_bm), ncol = ncol(raw_inp)-ncol(decontammed_data_low_bm) ) ) %>%
  diversity(index=div_ind)


diver_microdec <- rbind(microdec_df) %>% cbind( matrix(0,nrow=nrow(microdec_df), ncol = ncol(raw_inp)-ncol(microdec_df) ) ) %>%
  diversity(index=div_ind)

# decontammed_data %>%
# diversity(index=div_ind)


n_melanoma_samples <- length(diver_raw)
divers_df <- data.frame( 
  c( diver_raw, 
     diver_scrub, 
     diver_dec, diver_lb_dec, diver_rest, diver_microdec
  ), 
  c( rep('Raw', n_plasma_samples),
     rep('SCRuB', n_plasma_samples),
     rep('Decontam', n_plasma_samples),
     rep('Decontam (LB)', n_plasma_samples),
     rep('Restrictive', n_plasma_samples),
     rep('microDecon', n_plasma_samples) 
  )
)

colnames(divers_df) <- c('Shannon', 'Dataset')

divers_df %>% write.csv('../results/data/Fig3_plasma/plasma_shannon_diversities.csv')

set.seed(seed)


## following pipeline implemented by Poore et al

xy <-  metadataPSMatchedDPQCFiltered %>% 
  count(disease_type_consol) %>% 
  arrange(n) %>%
  pull(disease_type_consol) %>% as.character() %>%
  combn(2) 

#copy the sampling size's from the knight lab's notebook
# (https://github.com/biocore/tcga/blob/master/jupyter_notebooks/Plasma%20Kraken%20Machine%20Learning%20LOO%20Analysis.ipynb)
sampling_sizes <- c(25, 59, 69, 59, 69, 69)

dir.create(paste0('../results/data/Fig3_plasma/trial_', seed) )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/scrubbed_preds') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/raw_preds') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_preds') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_standard') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_lb') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/restrictive_preds') )
dir.create(paste0('../results/data/Fig3_plasma/trial_', seed, '/microdecon_preds') )

## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  scrubbed_hVsC_tmp <- loocvDTs(snmData = scrubbed_normalized,
                                samplingSize = sampling_sizes[idx], 
                                DTs = xy[ ,idx],
                                filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/scrubbed_preds/', xy[,idx][1], '_', xy[,idx][2] ),
                                caretTuneGrid = defaultGBMGrid)
}

## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  raw_hVsC_tmp <- loocvDTs(snmData = raw_normalized,
                           samplingSize = sampling_sizes[idx], 
                           DTs = xy[ ,idx],
                           filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/raw_preds/', xy[,idx][1], '_', xy[,idx][2] ),
                           caretTuneGrid = defaultGBMGrid)
}


## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  raw_hVsC_tmp <- loocvDTs(snmData = dec_normalized,
                           samplingSize = sampling_sizes[idx],
                           DTs = xy[ ,idx],
                           filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_preds/', xy[,idx][1], '_', xy[,idx][2] ),
                           caretTuneGrid = defaultGBMGrid)
}


## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  raw_hVsC_tmp <- loocvDTs(snmData = dec_standard_normalized,
                           samplingSize = sampling_sizes[idx], 
                           DTs = xy[ ,idx],
                           filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_standard/', xy[,idx][1], '_', xy[,idx][2] ),
                           caretTuneGrid = defaultGBMGrid)
}

## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  raw_hVsC_tmp <- loocvDTs(snmData = dec_lb_normalized, 
                           samplingSize = sampling_sizes[idx], 
                           DTs = xy[ ,idx],
                           filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/dec_lb/', xy[,idx][1], '_', xy[,idx][2] ),
                           caretTuneGrid = defaultGBMGrid)
}



## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  raw_hVsC_tmp <- loocvDTs(snmData = restrictive_normalized, 
                           samplingSize = sampling_sizes[idx], 
                           DTs = xy[ ,idx],
                           filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/restrictive_preds/', xy[,idx][1], '_', xy[,idx][2] ),
                           caretTuneGrid = defaultGBMGrid)
}



## loop through the 6 idifferent prediction tasks
for(idx in 1:ncol(xy)){
  scrubbed_hVsC_tmp <- loocvDTs(snmData = microdecon_normalized,
                                samplingSize = sampling_sizes[idx], 
                                DTs = xy[ ,idx],
                                filenameString = paste0('../results/data/Fig3_plasma/trial_', seed, '/microdecon_preds/', xy[,idx][1], '_', xy[,idx][2] ),
                                caretTuneGrid = defaultGBMGrid)
}



}









