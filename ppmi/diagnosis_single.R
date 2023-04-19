

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
mofa_imputed<-impute(MOFAobject, views=3)
mofa_imputed@imputed_data


group=1
vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]


### run the venn comparisons to load files 
# if we obtain the first 30 it means we are getting the top abslog2pval
# double check they are ordered like this!! 
n_top<-200
n_top_mir<-dim(signif_mirs)[1]
n_top_mir<-60
signif_rna_ord<-signif_rna[order(signif_rna[, 'padj']),]
dim(signif_rna)
dim(signif_mirs)[1]
signif_mirs<-signif_mirs[!is.na(signif_mirs$significant),]
signif_mirs
signif_mirs_ord<-signif_mirs[order(signif_mirs[, 'padj']),]
SIG_GENES<-signif_rna_ord[1:n_top,]$X
SIG_MIRNAS<-signif_mirs_ord[1:n_top_mir]$X
SIG_MIRNAS
SIG_PROTEINS<-signif_proteins$X



### then do predictions 
outdir_predict=paste0(outdir, '/predict/')
dir.create(outdir_predict)
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
yvar='COHORT'
yvar='COHORT'



choose_factors<-which(cors[,yvar]>0)
choose_factors
QUANTILE_THRESH<-0.02


detach(package:MOFA2,unload=TRUE)# replaces predict so we detach it 
require(MOFA2)
vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]

choose_factor_var<-vars_by_factor[choose_factors,]>1 
choose_factor_var
# choose factors based on important/associate

ws_all_miRNA<-get_weights(MOFAobject, factors=which(choose_factor_var[,1]), views=1)$miRNA
ws_all_RNA<-get_weights(MOFAobject, factors=which(choose_factor_var[,2]), views=2)$RNA
ws_all_prot<-get_weights(MOFAobject, factors=which(choose_factor_var[,3]), views=3)$proteomics


T1<-quantile(unlist(ws_all_RNA),QUANTILE_THRESH )
T2<-quantile(unlist(ws_all_miRNA), QUANTILE_THRESH)
T3<-quantile(unlist(ws_all_prot), QUANTILE_THRESH)
hist(ws_all_miRNA)
hist(ws_all_RNA)
hist(ws_all_prot)





### weight all 
### probably not a good idea to scale there might be high variability but NOT associated with disease control
# maybe better to scale by disease control association? 
#var_weights_miRNA<-ws_all_miRNA * vars_by_factor[choose_factors,'miRNA' ] /100
#var_weights_RNA<-ws_all_RNA * vars_by_factor[choose_factors,'RNA' ] / 100
all_feats_RNA_mofa<-rownames(ws_all_RNA)[rowSums((abs(ws_all_RNA)>abs(T1)))>0L]
all_feats_miRNA_mofa<-rownames(ws_all_miRNA)[rowSums((abs(ws_all_miRNA)>abs(T2)))>0L]
all_feats_prot_mofa<-rownames(ws_all_prot)[rowSums((abs(ws_all_prot)>abs(T3)))>0L]


length(all_feats_RNA_mofa)
length(all_feats_miRNA)


#hist(var_weights_miRNA[,1])
#hist(var_weights_miRNA[,2])
#get_data(MOFAobject, views='miRNA')[[1]]$group1
# these are the real datasets!!! 


######### Assign test data 
### TEST
test_data_miRNA<-assays(mofa_multi_complete_test)$miRNA
test_data_RNA<-assays(mofa_multi_complete_test)$RNA
test_data_prot<-assays(mofa_multi_complete_test)$proteomics
dim(test_data_miRNA); dim(test_data_RNA); 



y_actual<-colData(mofa_multi_complete_test)[,yvar]
length(y_actual)
CONF<-as.data.frame(colData(mofa_multi_complete_test)@listData[c('AGE_SCALED', 'SEX')])
CONF$SEX<-as.factor(CONF$SEX)

CONF_R<-t(CONF)




miRNA_data<-get_data(MOFAobject)$miRNA[[1]]
RNA_data<-get_data(MOFAobject)$RNA[[1]]
prot_data<-get_data(MOFAobject)$proteomics[[1]]
prot_data<-mofa_imputed@imputed_data$proteomics[[1]]


CONF_train<-MOFAobject@samples_metadata[c('AGE_SCALED', 'SEX')]
CONF_train$SEX<-as.factor(CONF_train$SEX)
CONF_train_R<-t(CONF_train)


#######################################################
# here using DE GENES 

### FEATURE SELECTION

#### ADD ALSO SEX AND AGE 
use_mofa=FALSE
if (use_mofa){
  
  RNA_data
  all_feats_RNA=all_feats_RNA_mofa
  all_feats_miRNA=all_feats_miRNA_mofa
  all_feats_prot=all_feats_prot_mofa
}else{
  all_feats_RNA=SIG_GENES
  all_feats_miRNA=SIG_MIRNAS
  all_feats_prot=SIG_PROTEINS
}


###### filter by significant genes/ features 
test_data_miRNA_filt<-test_data_miRNA[ rownames(test_data_miRNA) %in% all_feats_miRNA,]
test_data_RNA_filt<-test_data_RNA[ rownames(test_data_RNA) %in% all_feats_RNA,]
test_data_prot_filt<-test_data_prot[ rownames(test_data_prot) %in% all_feats_prot,]

length(SIG_PROTEINS_FILT);length(all_feats_RNA)
length(SIG_miRNA_FILT);length(all_feats_miRNA)


dim(test_data_miRNA_filt);dim(test_data_RNA_filt);dim(test_data_prot_filt)


# Filter by a smaller set that does not contain NA OR impute!!!
remove_na<-TRUE
if (remove_na){
    # FILTER VALS WITH NAS 
    df=test_data_prot_filt; SIG_PROTEINS_FILT<-names(which(rowSums(is.na(df)) ==0))
    test_data_prot_filt<-test_data_prot[ rownames(test_data_prot) %in% SIG_PROTEINS_FILT,]
    
    df=test_data_RNA_filt; SIG_RNA_FILT<-names(which(rowSums(is.na(df)) ==0));
    test_data_RNA_filt<-test_data_RNA[ rownames(test_data_RNA) %in% SIG_RNA_FILT,]
    
    df=test_data_miRNA_filt; SIG_miRNA_FILT<-names(which(rowSums(is.na(df)) ==0));
    test_data_miRNA_filt<-test_data_miRNA[ rownames(test_data_miRNA) %in% SIG_miRNA_FILT,]
  
  
  
  prot_data_filt<-prot_data[rownames(prot_data) %in% unique(SIG_PROTEINS_FILT),]
  RNA_data_filt<-RNA_data[rownames(RNA_data) %in% unique(SIG_RNA_FILT),]; dim(RNA_data_filt)
  miRNA_data_filt<-miRNA_data[rownames(miRNA_data) %in% unique(SIG_miRNA_FILT),]
  
  
}else{
  RNA_data_filt<-RNA_data[rownames(RNA_data) %in% unique(all_feats_RNA),]; dim(RNA_data_filt)
  miRNA_data_filt<-miRNA_data[rownames(miRNA_data) %in% unique(all_feats_miRNA),]
  prot_data_filt<-prot_data[rownames(prot_data) %in% unique(all_feats_prot),]
  
}


dim(RNA_data_filt);dim(miRNA_data_filt);dim(prot_data_filt)


dim(RNA_data)
dim(RNA_data_filt);  dim(test_data_RNA_filt); 
dim(miRNA_data_filt) ;dim(test_data_miRNA_filt); 

run_single=TRUE;
single_mode='protein';
single_mode='miRNA';
single_mode='RNA'

if (run_single){
    if (single_mode=='RNA'){
      data_filt<-RNA_data_filt
      test_data_filt<-test_data_RNA_filt
    }else if (single_mode=='miRNA'){
      data_filt<-miRNA_data_filt
      test_data_filt<-test_data_miRNA_filt
    }else if (single_mode=='protein'){
      data_filt<-prot_data_filt
      test_data_filt<-test_data_prot_filt
    }
  
}else{
  print('run ALL modalities')
  data_filt<-rbind(RNA_data_filt, miRNA_data_filt, prot_data_filt)
  test_data_filt<-rbind(test_data_RNA_filt, test_data_miRNA_filt, test_data_prot_filt)
}




rownames(data_filt); rownames(test_data_filt)

## append confounding factors 

test_data_filt_conf<-rbind(test_data_filt, CONF_R)
data_filt_conf<-rbind(data_filt, CONF_train_R)
test_data_filt_conf['SEX',]
data_filt_conf['SEX',]


dim(data_filt)

rlang::env_unlock(env = asNamespace('caret'))

library(randomForest)

detach(package:MOFA2,unload=TRUE)# replaces predict so we detach it 


#if (add_conf){
#  all_preds(data_filt_conf, test_data_filt_conf)
  
#}else{
#  stats<-all_preds(data_filt, test_data_filt)
#}


add_conf=FALSE

#all_preds<-function(data_filt, test_data_filt){
if (add_conf){
  data_filt_in=data_filt_conf
  test_data_filt_in=test_data_filt_conf
}else{
  data_filt_in=data_filt
  test_data_filt_in=test_data_filt
}
rownames(test_data_filt)


add_conf
  # Prepare data
  # Predict EORTC.risk with factor 1,2 only!
  df <- as.data.frame(t(data_filt_in))
  colnames(df)<-gsub('-', '.',colnames(df) ) 
  
  rownames(test_data_filt_in)<-gsub('-','.', rownames(test_data_filt_in))
  df_test <- as.data.frame(t(test_data_filt_in))
  
  
  colnames(df_test)==  colnames(df)
  dim(df_test)
  # Train the model for IGHV
  
  y_predict='CONCOHORT_DEFINITION'
  y_predict=yvar
 # as.factor(MOFAobject@samples_metadata[,'CONCOHORT_DEFINITION'])
  
  
  y <- as.factor(MOFAobject@samples_metadata[,y_predict])
  
  y
  #run_rf( df, ){
  #model.y <- randomForest(y ~ .,data= df, ntree=35)
  #tuning_res<-tuneRF(df, y, ntreeTry = 500, doBest = TRUE, stepFactor = 2  )
  tune_res=TRUE
  if (tune_res){
    x=df
    df;y
    ntry=10
    repeats=3
    control <- trainControl(method="repeatedcv", number=ntry, repeats=repeats)
    
    seed <- 30
    metric <- "Accuracy"
    set.seed(seed)
    mtry <- floor(sqrt(ncol(x))) ### mtry will depend on the number of features 
    tunegrid <- expand.grid(.mtry=mtry); 
    print(paste0('gridsize: ',as.character(tunegrid$.mtry)))
    df_all<-df;df_all$y=y
    method='cforest'## cforest takes longer, but it looks like it has better outcome
    method='rf'## cforest takes longer, but it looks like it has better outcome
    
    
    rf_default <- train(y~., data=df_all, method='cforest', metric=metric, 
                        tuneGrid=tunegrid, trControl=control)
    rf_default
    
    
  
    print(rf_default)
  }
  #model.y_tuned<-tuning_res
  tune_set='default'
  
  if (tune_set=='random'){
    #rf.imputed <- rfImpute(x = x, y=y)
    
    rf_random <- train(x = x, y=y, method=method, metric=metric,
                       tuneLength=10, trControl=control)
    model.y_tuned<-rf_random
    
  }else{
    model.y_tuned<-rf_default
    
  }
  
  model.y_tuned$finalModel
 
  # Do predictions
    
    ###
  #head(test_data_filt_in)
  #y.pred <- stats::predict(model.y, df_test)
  unique(rownames(df_test))
  colnames(df_test)
  y.pred <- stats::predict(model.y_tuned, df_test)
  length(y.pred)
  dim(df_test)
  df_test
  colnames(df)
  #### ON TEST SET 
  colnames(df_test)
  colnames(df)
  predicted <- y.pred
  actual <-as.factor(y_actual)
  predicted; actual
  
  confusion_mat = as.matrix(table(actual, predicted )) 
  #predictions<-as.data.frame(cbind(c(actual), c(predicted)))
  #colnames(predictions)=c('observed', 'predicted')
  conf_mat<-confusionMatrix(confusion_mat)
  print(confusion_mat)
  conf_mat$byClass
  print(conf_mat$byClass['Balanced Accuracy'])
  
  
  important_feats<-round(importance(model.y_tuned$finalModel@weights), 2)
  important_feats<-order(-important_feats)
  params_train<-paste0(outdir_predict, 'predict_', yvar, use_mofa, '_t_', QUANTILE_THRESH, 
                 'add', 'n_top', n_top,'n_top_mir', n_top_mir, add_conf , 'single_', run_single ,
                 'smode_', single_mode,  '_tune_', ntry, repeats, 
                 mtry, tune_set, seed, seed_tr_test)
  write.csv(important_feats,  paste0( params_train, '_important.csv'))
  write.csv(conf_mat$byClass, paste0(params_train,   '_results.csv'))

  length(dim(conf_mat$byClass))
  if (length(dim(conf_mat$byClass))){
    accuracy_w=conf_mat$byClass[,'Balanced Accuracy']; 
    accuracy_w=mean(accuracy_w, na.rm=TRUE)
  }else{
    accuracy_w=conf_mat$byClass['Balanced Accuracy']; 
    
  }
  df_stats=  c(yvar, use_mofa, QUANTILE_THRESH,n_top, n_top_mir,
                     add_conf, run_single,single_mode, method, ntry , repeats, mtry,tune_set,
                     seed, seed_tr_test, accuracy_w)
  
  
  write.table(t(df_stats), paste0(outdir,'all_stats.csv'), append=TRUE,sep=',', col.names = FALSE)
  print(accuracy_w)
  #return(conf)
#}






#library(pROC)
#library(pROC)
#
#roc.mock <- roc(ifelse(predictions$observed==3, 3, 2), as.numeric(predictions$predicted))
#plot(roc.mock, col = "gray60")
#


  
  
  # example 2-class data and model training
#install.packages('caret')
#library(caret)
#d <- iris[51:150,]
#d[,5] <- factor(d[,5])
#model <- train(x = d[,c(1,3)], y = d[,5], method = 'lda', metric = 'ROC', 
#               trControl=trainControl(method = 'repeatedcv', number = 10, 
#                                      repeats = 10, savePredictions = T, 
#                                      classProbs = T, summary = twoClassSummary))
#
#

print(confusion_mat)
  
print(conf_mat$byClass['Balanced Accuracy'])
  
