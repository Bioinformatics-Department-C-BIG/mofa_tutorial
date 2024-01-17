suppressPackageStartupMessages(library(randomForest))
packages <- c('dplyr', 'ggplot2', 'caret')
invisible(lapply(packages, library, character.only = TRUE))

suppressPackageStartupMessages(library(randomForest))







run_train_validation<-function(df){
  ## Define repeated cross validation with 5 folds and three repeats
  repeat_cv <- trainControl(method='repeatedcv', number=5, repeats=3)
  train_index <- createDataPartition(y=df$y, p=0.7, list=FALSE)
  val_folds<-createFolds(df$y, k = 10, list = TRUE, returnTrain = TRUE)
  
  res<-sapply(val_folds, cross_val_score, df=df)
  return(res)
  
}


cross_val_score<- function(train_index, df){
  #'
  #' @train_index 
  #'
  #'
  
  repeat_cv <- trainControl(method='repeatedcv', number=5, repeats=3)
  
  training_set <- df[train_index, ]
  testing_set <- df[-train_index, ]
  
  
  
  
  mtry <- floor(sqrt(ncol(df)+1)) ### mtry will depend on the number of features 
  tunegrid <- expand.grid(.mtry=mtry);
  
  
  ##### tuning on training 
  rf_default <- train(y~., data=df, method='cforest', 
                      trControl=repeat_cv,tuneGrid=tunegrid )
  rf_default
  # Do predictions
  
  testing_set$pred <- stats::predict(rf_default, testing_set)
  
  predicted<-as.factor(testing_set$pred)
  actual = testing_set$y
  confusion_mat = as.matrix(table(actual, predicted )) 
  confusion_mat
  confusion_mat
  acc<-caret::confusionMatrix(as.matrix(table(actual, predicted )) )
  return(acc$byClass)
}
