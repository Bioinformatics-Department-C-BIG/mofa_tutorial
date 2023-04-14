

library(caret)
#remove.packages('rlang')
#install.packages('vctrs')
library('vctrs')

library(randomForest)
#install.packages('mlbench')
library(mlbench)
library(caret)


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
SIG_GENES
### then do predictions 
T=0.3;

############### MOFA SPECIFIC WEIGHTED GENES ##################################
###############################################################################
# Should I overlap with de? 
# choose factors with higher overlap with clinvar of interest 
###############################################################################

cors<-correlate_factors_with_covariates(MOFAobject,
                                        covariates = names(non_na_vars), 
                                        plot = "log_pval", 
                                        return_data = TRUE
                                        
)
yvar='CONCOHORT'
choose_factors<-which(cors[,yvar]>0)
choose_factors
detach(package:MOFA2,unload=TRUE)# replaces predict so we detach it 
require(MOFA2)
# choose factors based on important/associate
ws_all_miRNA<-get_weights(MOFAobject, factors=choose_factors)$miRNA 
ws_all_RNA<-get_weights(MOFAobject, factors=choose_factors)$RNA 

QUANTILE_THRESH<-0.05
T1<-quantile(ws_all_RNA,QUANTILE_THRESH )
T2<-quantile(ws_all_miRNA, QUANTILE_THRESH)
hist(ws_all_miRNA)






### weight all 
### probably not a good idea to scale there might be high variability but NOT associated with disease control
# maybe better to scale by disease control association? 
#var_weights_miRNA<-ws_all_miRNA * vars_by_factor[choose_factors,'miRNA' ] /100
#var_weights_RNA<-ws_all_RNA * vars_by_factor[choose_factors,'RNA' ] / 100
all_feats_RNA<-rownames(ws_all_RNA)[rowSums((abs(ws_all_RNA)>abs(T1)))>0L]
all_feats_miRNA<-rownames(ws_all_miRNA)[rowSums((abs(ws_all_miRNA)>abs(T2)))>0L]

length(all_feats_RNA)
length(all_feats_miRNA)


#hist(var_weights_miRNA[,1])
#hist(var_weights_miRNA[,2])
#get_data(MOFAobject, views='miRNA')[[1]]$group1
# these are the real datasets!!! 



miRNA_data<-get_data(MOFAobject)$miRNA[[1]]
RNA_data<-get_data(MOFAobject)$RNA[[1]]
CONF_train<-MOFAobject@samples_metadata[c('AGE_AT_VISIT', 'SEX')]
CONF_train$SEX<-as.factor(CONF_train$SEX)
CONF_train_R<-t(CONF_train)



### TEST
test_data_miRNA<-assays(mofa_multi_complete_test)$miRNA
test_data_RNA<-assays(mofa_multi_complete_test)$RNA
y_actual<-colData(mofa_multi_complete_test)[,yvar]
CONF<-as.data.frame(colData(mofa_multi_complete_test)@listData[c('AGE_AT_VISIT', 'SEX')])
CONF$SEX<-as.factor(CONF$SEX)

CONF_R<-t(CONF)
#######################################################
# here using DE GENES 

### FEATURE SELECTION

#### ADD ALSO SEX AND AGE 
use_mofa=TRUE
if (use_mofa){
  
  all_feats_RNA
  RNA_data
  
}else{
  all_feats_RNA=SIG_GENES
}




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
test_data_filt_conf['SEX',]
data_filt_conf['SEX',]


dim(data_filt)

rlang::env_unlock(env = asNamespace('caret'))

library(randomForest)

detach(package:MOFA2,unload=TRUE)# replaces predict so we detach it 


all_preds<-function(data_filt, test_data_filt){

  # Prepare data
  # Predict EORTC.risk with factor 1,2 only!
  df <- as.data.frame(t(data_filt))
  colnames(df)<-gsub('-', '.',colnames(df) ) 
  
  
  rownames(test_data_filt)<-gsub('-','.', rownames(test_data_filt))
  df_test <- as.data.frame(t(test_data_filt))
  colnames(df_test)==  colnames(df)
  # Train the model for IGHV
  
  y_predict='CONCOHORT_DEFINITION'
  y_predict=yvar
 # as.factor(MOFAobject@samples_metadata[,'CONCOHORT_DEFINITION'])
  
  
  y <- as.factor(MOFAobject@samples_metadata[,y_predict])
  
  y
  #run_rf( df, ){
  model.y <- randomForest(y ~ .,data= df, ntree=35)
  #tuning_res<-tuneRF(df, y, ntreeTry = 500, doBest = TRUE, stepFactor = 2  )
  tune_res=TRUE
  if (tune_res){
    x=df
    df;y
    ntry=5
    repeats=2
    control <- trainControl(method="repeatedcv", number=ntry, repeats=repeats)
    
    seed <- 7
    metric <- "Accuracy"
    set.seed(seed)
    mtry <- floor(sqrt(ncol(x))/5)
    
    tunegrid <- expand.grid(.mtry=mtry)
    df_all<-df;df_all$y=y
    #rf_default <- train(y~., data=df_all, method="rf", metric=metric, 
    #                    tuneGrid=tunegrid, trControl=control)
    rf_random <- caret::train(x = x, y=y, method="cforest", metric=metric,
                              tuneLength=15, trControl=control)
    
  
    print(rf_default)
  }
  #model.y_tuned<-tuning_res
  model.y_tuned<-rf_random

 
  # Do predictions
    
    ###
  head(test_data_filt)
  #y.pred <- stats::predict(model.y, df_test)
  
  y.pred <- stats::predict(model.y_tuned, df_test)

  #### ON TEST SET 
  colnames(df_test)
  colnames(df)
  predicted <- y.pred
  actual <-as.factor(y_actual)
  confusion_mat = as.matrix(table(actual, predicted )) 
  #predictions<-as.data.frame(cbind(c(actual), c(predicted)))
  #colnames(predictions)=c('observed', 'predicted')
  conf<-confusionMatrix(confusion_mat)
  print(conf$byClass)
  
  print(confusion_mat)
  round(importance(model.y), 2)
  
  write.csv(conf$byClass, paste0(outdir, 'predict_', yvar, use_mofa, '_t_', QUANTILE_THRESH, 'add', add_conf  ,  '.csv'))
  
  
  return(confusion_mat)
}
add_conf=TRUE
if add_conf{
  all_preds(data_filt_conf, test_data_filt_conf)
  
}else{
  all_preds(data_filt, test_data_filt)
}





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






