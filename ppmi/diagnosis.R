
fetch_top_weights<-function(ws){
  ws=ws[order(-ws)]
  ws_top<-ws[abs(ws)>T]
  #return(ws_top)
  return(names(ws_top))
}

T=0.5; ## WEIGHT BY VARIANCE CAPTURED, RANK IN FACTOR, WEIGHT ETC. 
ws_all_miRNA<-get_weights(MOFAobject, factors=c(2,3,4,6))$miRNA 
ws_all_RNA<-get_weights(MOFAobject, factors=c(2,3,4,6))$RNA 

all_feats_miRNA<-unlist(apply(ws_all_miRNA,2,fetch_top_weights) )
all_feats_RNA<-unlist(apply(ws_all_RNA,2,fetch_top_weights) )
all_feats_miRNA

miRNA_data<-get_data(MOFAobject)$miRNA[[1]]
RNA_data<-get_data(MOFAobject)$RNA[[1]]


### FEATURE SELECTION

RNA_data_filt<-RNA_data[unique(all_feats_RNA),]
miRNA_data_filt<-miRNA_data[unique(all_feats_miRNA),]

dim(RNA_data_filt)
dim(miRNA_data_filt)


data_filt<-rbind(RNA_data_filt, miRNA_data_filt)
dim(data_filt)




# Prepare data
# Predict EORTC.risk with factor 1,2 only!
df <- as.data.frame(t(data_filt))
# Train the model for IGHV
y_predict='CONCOHORT_DEFINITION'
colnames(df)<-gsub('-', '.',colnames(df) ) 
df$y <- as.factor(MOFAobject@samples_metadata[,y_predict])
model.y <- randomForest(y ~ .,data= df, ntree=10)
df$y <- NULL # important 


# Do predictions
MOFAobject@samples_metadata$y.pred <- stats::predict(model.y, df)
MOFAobject@samples_metadata$y.pred

# Assess performance 
predicted<-MOFAobject@samples_metadata$y.pred
actual <-as.factor(MOFAobject@samples_metadata[,y_predict])
confusion_mat = as.matrix(table(actual, predicted )) 

print(confusion_mat)
View(confusion_mat)
round(importance(model.y), 2)



