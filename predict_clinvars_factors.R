
######### 1. Predict cohort variables ####

######### 1. Predict clinical variables ####
# Input mofa or nmf factors in it
#### and label ####

suppressPackageStartupMessages(library(randomForest))
packages <- c('dplyr', 'ggplot2', 'caret', 'party')
invisible(lapply(packages, library, character.only = TRUE))

suppressPackageStartupMessages(library(randomForest))



source(paste0(script_dir, 'ppmi/predict_utils.R'))


# Prepare data
df_mofa <- as.data.frame(get_factors(MOFAobject, factors=sel_factors)[[1]])

df_mofa <- as.data.frame(get_factors(MOFAobject, factors=1:15)[[1]])







# Do predictions with the factors 
# Train the model for eortc.risk

# Use the pricnipal components 
# Train-validation split #
# Cross-Validation ####


df_mofa$y<- as.factor(MOFAobject@samples_metadata$COHORT)
df_mofa_age <- cbind(df_mofa,MOFAobject@samples_metadata[, c('AGE_SCALED', 'SEX')])
res_age_mofa<-run_train_validation( df=df_mofa_age)



#### df2: NMF ####
##################
rownames(h_t);
df1$PATNO
df_nmf <- h_t
df_nmf$y = as.factor(df1$COHORT)
df_age_nmf <- cbind(df_nmf,df1[, c('AGE_SCALED', 'SEX')])



## Set seed for reproducibility
set.seed(123)






## TODO: issues : there is class imbalance 
# TODO: issues 

res_age_nmf<-run_train_validation( df=df_age_nmf)
res_age_mofa
res<-run_train_validation(df_nmf)
res_age<-run_train_validation(df=df_age)



mean(res[ 'Balanced Accuracy',])  
mean(res_age[ 'Balanced Accuracy',])  
mean(res_age_nmf[ 'Balanced Accuracy',], na.rm=TRUE)  ; sd(res_age_nmf[ 'Balanced Accuracy',], na.rm = TRUE)  
mean(res_age_mofa[ 'Balanced Accuracy', na.rm=TRUE]) ;sd(res_age_mofa[ 'Balanced Accuracy',])  
 

res  



