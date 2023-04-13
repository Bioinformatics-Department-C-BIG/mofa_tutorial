

library(caret)
#remove.packages('rlang')
#install.packages('vctrs')
library('vctrs')
library(caret)
install.packages('ConfusionTableR')
library(ConfusionTableR)

#update.packages("rlang")

suppressPackageStartupMessages(library(randomForest))

#### make sure that mofa was ran with complete cases otherwise there might be issues 

#install.packages('pROC')

#### FIRST obtain the highly weighted features from MOFA 

fetch_top_weights<-function(ws){
  #ws=ws_all_miRNA
  ### apply to each factor 
  #ws=ws_all_miRNA[,1]
  ws=ws[order(-ws)]
  ws_top<-ws[abs(ws)>T]
  ws_top
  #return(ws_top)
  return(names(ws_top))
}


sig<-read.csv(paste0(outdir_s, '/significant0.005_0.25.csv'))
SIG_GENES<-sig$X

### then do predictions 
T=0.5;
all_preds<-function(T){
  ## WEIGHT BY VARIANCE CAPTURED, RANK IN FACTOR, WEIGHT ETC. 
  ## here automatically obtain highly associated features 
  ws_all_miRNA<-get_weights(MOFAobject, factors=c(2,3,4,6))$miRNA 
  ws_all_RNA<-get_weights(MOFAobject, factors=c(2,3,4,6))$RNA 
  
  ws_all_miRNA * vars_by_factor
  
  
  all_feats_miRNA<-unlist(apply(ws_all_miRNA,2,fetch_top_weights) )
  all_feats_RNA<-unlist(apply(ws_all_RNA,2,fetch_top_weights) )
  all_feats_miRNA
  
  miRNA_data<-get_data(MOFAobject)$miRNA[[1]]
  RNA_data<-get_data(MOFAobject)$RNA[[1]]
  
  
  ### FEATURE SELECTION
  sel_feats<-unique(all_feats_RNA)
  
  select_features<-function(select_feats){
    RNA_data_filt<-RNA_data[SIG_GENES,]
    
    
  }
  
  
 # miRNA_data_filt<-miRNA_data[unique(all_feats_miRNA),]
  
#  dim(RNA_data_filt)
  #dim(miRNA_data_filt)
  
  
 # data_filt<-rbind(RNA_data_filt, miRNA_data_filt)
  data_filt<-RNA_data_filt
  dim(data_filt)
  
  
  
  
  # Prepare data
  # Predict EORTC.risk with factor 1,2 only!
  df <- as.data.frame(t(data_filt))
  # Train the model for IGHV
  y_predict='CONCOHORT_DEFINITION'
  y_predict='CONCOHORT'
  
  
  colnames(df)<-gsub('-', '.',colnames(df) ) 
  
  df$y <- as.factor(MOFAobject@samples_metadata[,y_predict])
  
  
  #run_rf( df, ){
    model.y <- randomForest(y ~ .,data= df, ntree=10)
    df$y <- NULL # important 
    
  #}
 
  
  
  
  
  # Do predictions
  MOFAobject@samples_metadata$y.pred <- stats::predict(model.y, df)
  MOFAobject@samples_metadata$y.pred
  
  # Assess performance 
  ## diagnostic
 # install.packages('caret')
  
  predicted<-MOFAobject@samples_metadata$y.pred
  actual <-as.factor(MOFAobject@samples_metadata[,y_predict])
  confusion_mat = as.matrix(table(actual, predicted )) 
  predictions<-as.data.frame(cbind(c(actual), c(predicted)))
  colnames(predictions)=c('observed', 'predicted')
  #conf<-confusionMatrix(confusion_mat)
  
  accuracy(confusion_mat)
  print(confusion_mat)
  round(importance(model.y), 2)
  
  mc_df <- ConfusionTableR::multi_class_cm(predictions$actual, 
                                           predictions$predicted,
                                           mode="everything")
  
  return(conf)
}

all_preds(T=0.1)

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






