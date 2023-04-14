

library(caret)
#remove.packages('rlang')
#install.packages('vctrs')
library('vctrs')

library(caret)
#install.packages('ConfusionTableR')
library(ConfusionTableR)

sessionInfo()
library(dplyr)
library(randomForest) 
#install.packages('Metrics')
library(ranger)
library(Metrics)


#update.packages("rlang")

suppressPackageStartupMessages(library(randomForest))

#### make sure that mofa was ran with complete cases otherwise there might be issues 

#install.packages('pROC')

#### FIRST obtain the highly weighted features from MOFA 


vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]



sig<-read.csv(paste0(outdir_s, '/significant0.005_0.25.csv'))
SIG_GENES<-sig$X

### then do predictions 
T=0.3;

############### MOFA SPECIFIC WEIGHTED GENES ##################################
###############################################################################
# Should I overlap with de? 
# choose factors with higher overlap with clinvar of interest 
###############################################################################

choose_factors<-c(2,3)
# choose factors based on important/associate
ws_all_miRNA<-get_weights(MOFAobject, factors=choose_factors)$miRNA 
ws_all_RNA<-get_weights(MOFAobject, factors=choose_factors)$RNA 

QUANTILE_THRESH<-0.5
T1<-quantile(ws_all_RNA,QUANTILE_THRESH )
T2<-quantile(ws_all_miRNA, QUANTILE_THRESH)
hist(ws_all_miRNA)






### weight all 
### probably not a good idea to scale there might be high variability but NOT associated with disease control
# maybe better to scale by disease control association? 
#var_weights_miRNA<-ws_all_miRNA * vars_by_factor[choose_factors,'miRNA' ] /100
#var_weights_RNA<-ws_all_RNA * vars_by_factor[choose_factors,'RNA' ] / 100
head(var_weights_RNA)
all_feats_RNA<-rownames(ws_all_RNA)[rowSums((abs(ws_all_RNA)>abs(T1)))>0L]
all_feats_miRNA<-rownames(ws_all_miRNA)[rowSums((abs(ws_all_miRNA)>abs(T2)))>0L]

length(all_feats_RNA)
length(all_feats_miRNA)


summary(var_weights_miRNA)
#hist(var_weights_miRNA[,1])
#hist(var_weights_miRNA[,2])
#get_data(MOFAobject, views='miRNA')[[1]]$group1
# these are the real datasets!!! 



miRNA_data<-get_data(MOFAobject)$miRNA[[1]]
RNA_data<-get_data(MOFAobject)$RNA[[1]]
CONF_train<-MOFAobject@samples_metadata[c('AGE_AT_VISIT', 'SEX')]
CONF_train_R<-t(CONF_train)



### TEST
test_data_miRNA<-assays(mofa_multi_complete_test)$miRNA
test_data_RNA<-assays(mofa_multi_complete_test)$RNA
y_actual<-colData(mofa_multi_complete_test)[,'CONCOHORT']
CONF<-as.data.frame(colData(mofa_multi_complete_test)@listData[c('AGE_AT_VISIT', 'SEX')])
CONF_R<-t(CONF)
#######################################################
# here using DE GENES 

### FEATURE SELECTION

#### ADD ALSO SEX AND AGE 
use_mofa=TRUE
if (use_mofa){
  RNA_data
  RNA_data_filt<-RNA_data[rownames(RNA_data) %in% unique(all_feats_RNA),]
  miRNA_data_filt<-miRNA_data[rownames(miRNA_data) %in% unique(all_feats_miRNA),]
  
  
  
  
  test_data_miRNA_filt<-test_data_miRNA[ rownames(test_data_miRNA) %in% all_feats_miRNA,]
  test_data_RNA_filt<-test_data_RNA[ rownames(test_data_RNA) %in% all_feats_RNA,]
  
  
  
  dim(RNA_data)
  dim(RNA_data_filt);  dim(test_data_RNA_filt); 
  dim(miRNA_data_filt)
  
  
  data_filt<-rbind(RNA_data_filt, miRNA_data_filt)
  test_data_filt<-rbind(test_data_RNA_filt, test_data_miRNA_filt)
  
  ## append confounding factors 
  
  test_data_filt_conf<-rbind(test_data_filt, CONF_R)
  data_filt_conf<-rbind(data_filt, CONF_train_R)
  
  
  
  
  dim(data_filt)
}else{
  data_filt<-RNA_data_filt
  
}

#sel_feats<-unique(all_feats_RNA)


select_features<-function(select_feats){
  RNA_data_filt<-RNA_data[SIG_GENES,]
  
  
}


all_preds<-function(data_filt){
 
  # Prepare data
  # Predict EORTC.risk with factor 1,2 only!
  df <- as.data.frame(t(data_filt))
  colnames(df)<-gsub('-', '.',colnames(df) ) 
  
  
  rownames(test_data_filt)<-gsub('-','.', rownames(test_data_filt))
  df_test <- as.data.frame(t(test_data_filt))
  
  # Train the model for IGHV
  y_predict='CONCOHORT_DEFINITION'
  y_predict='CONCOHORT'
 # as.factor(MOFAobject@samples_metadata[,'CONCOHORT_DEFINITION'])
  
  
  df$y <- as.factor(MOFAobject@samples_metadata[,y_predict])
  
  
  #run_rf( df, ){
  model.y <- randomForest(y ~ .,data= df, ntree=35)
  df$y <- NULL # important 
    
  #}
 
  # Do predictions
    
    ###
  head(test_data_filt)
  
  y.pred <- stats::predict(model.y, df_test)

  
  
  #### ON TEST SET 
  
  # Assess performance 
  ## diagnostic
 # install.packages('caret')
  
  predicted <- y.pred
  actual <-as.factor(y_actual)
  confusion_mat = as.matrix(table(actual, predicted )) 
  #predictions<-as.data.frame(cbind(c(actual), c(predicted)))
  #colnames(predictions)=c('observed', 'predicted')
  #conf<-confusionMatrix(confusion_mat)

  
  print(confusion_mat)
  round(importance(model.y), 2)
  
  
  
  return(confusion_mat)
}

all_preds(data_filt_conf)



library(pROC)
library(pROC)

roc.mock <- roc(ifelse(predictions$observed==3, 3, 2), as.numeric(predictions$predicted))
plot(roc.mock, col = "gray60")


# example 2-class data and model training
#install.packages('caret')
library(caret)
d <- iris[51:150,]
d[,5] <- factor(d[,5])
model <- train(x = d[,c(1,3)], y = d[,5], method = 'lda', metric = 'ROC', 
               trControl=trainControl(method = 'repeatedcv', number = 10, 
                                      repeats = 10, savePredictions = T, 
                                      classProbs = T, summary = twoClassSummary))






