install.packages('e1071')
library(e1071)


y<-as.factor(Y_raw$Subtype)
  
  
  glm.fit(y, dat )
  
x=cbind(genes_sel_features, proteins_sel_features)
  
dat=data.frame( y=y,x=x)
  


train=dat[1:15,];y_train=y[1:15]
test=dat[16:16,];y_test=y[16:16]
svmfit = svm( ~ ., data = train, kernel = "linear", cost = 10, scale = FALSE)


y_test
predict(svmfit, test)
        
