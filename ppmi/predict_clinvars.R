
library(GGally)
source(paste0(script_dir, 'ppmi/utils.R'))

#### Make the selection elsewhere and input the features here
cors_both<-get_correlations(MOFAobject = MOFAobject, covariates = 'CONCOHORT')

sel_factors=which(cors_both[[1]]>-log10(0.05))

###

### Selected features  - taken from VISIT 8 MOFA for training
# TODO: use tuning to decide how many
# TODO: compare miRNAs and RNAs 
w_fs<-get_weights(MOFAobject, view='RNA', factors = c(4))[[1]]
selected_feats<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.995),])
x_rnas<-MOFAobject@data$RNA[[1]][selected_feats, ]

w_fs<-get_weights(MOFAobject, view='miRNA', factors = c(3))[[1]]
selected_feats2<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.95),])
w_fs<-get_weights(MOFAobject, view='miRNA', factors = c(4))[[1]]
selected_feats3<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.95),])
selected_feats2<-c(selected_feats2,selected_feats3 )

selected_feats2
length(selected_feats)
####


### SELECT FACTORS ####
Z <- get_factors(MOFAobject)[[1]]
x_factors<-data.frame(cbind(Z[,c(sel_factors)], y))
x_factors
GGally::ggpairs(x_factors)
x_vars<-cbind(Z[,c(sel_factors, 15)],MOFAobject@samples_metadata[, c('AGE', 'SEX')]   )
x_vars<-cbind(Z[,c(sel_factors, 15)])
x_vars<-cbind(Z[,c(sel_factors, 15)],MOFAobject@samples_metadata[, c('AGE')]   )
x_vars<-cbind(Z[,c(sel_factors, 15)]  )

x_vars


#### SELECT MOLECULES  from baseline objects ####
# 1. Load baseline objects - after VSN- se_filt?
process_mirnas=TRUE
source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
se_filt_V06<-filter_se(se, VISIT='V06', sel_coh,sel_ps)
se_filt_V04<-filter_se(se, VISIT='V04', sel_coh,sel_ps)

### Selected features from baseline
# TODO: split by caret package 
library(caret)
#install.packages('vctrs')
set.seed(34569)
# split so that there are equal numbers of different stages of NHY 
# For categorical use cut function
common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )

x_vars_sel=preprocess_visit_predict(se_filt_BL, common=common,feat_names = selected_feats2)
trainIndex <- createDataPartition(x_vars_sel$NHY, p = .8, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)  
x_vars_sel$CONCOHORT
hist(x_vars_sel$NP3_TOT)
cut(x_vars_sel$NP3_TOT, breaks = 4)


trainIndex=1:floor(dim(x_vars_sel)[1] *0.7)
trainIndex
x_vars_sel_train<- x_vars_sel[trainIndex, ]
x_vars_sel_test<- x_vars_sel[-trainIndex, ]

x_vars_sel_train$SEX
x_vars_sel_test$SEX

#######






#### Outcome is from Visit 8
#### 1. V
#
#
# TODO: train test split: 
# 1. subset by samples, and by features
# ALSO split to train test
library(glmnet)
y= se_filt_V08[,match(x_vars_sel_train$PATNO,se_filt_V08$PATNO )]$NP3_TOT
y_test<-se_filt_V08[,match(x_vars_sel_test$PATNO,se_filt_V08$PATNO )]$NP3_TOT
head(x_vars_sel)
dim(x_vars)
dim(x_vars)


## PREPROCESSING 


x_vars_sel_train

## Molecular markers convert to LOGFC
get_logFC_bySample<-function(x_vars,feature){
  #' for each sample
  cn_mean<-mean(x_vars_controls[,feature ])
  return((x_vars[,feature]-cn_mean)/cn_mean)
  
  
}
to_Add<-colnames(x_vars_sel_train[1:20])
x_vars_controls<-x_vars_sel_train[x_vars_sel_train$CONCOHORT==2,]

log_FCs<-as.data.frame(sapply(to_Add,get_logFC_bySample,x= x_vars_sel_train))

x_vars_sel_train_lfc<-mutate(x_vars_sel_train, log_FCs)

x_vars_sel_train_lfc<-x_vars_sel_train




library(dplyr)
library(tidyverse)


scale_numeric<-function(x_vars){
  #''DO not scale y? 
  #'
  x_scaled<-as.data.frame(as_tibble( x_vars) %>%
    dplyr::select(-c(PATNO_EVENT_ID, PATNO))%>%
    mutate_all(as.numeric)%>%
    
    # TODO: scale with caret/based on training
    mutate_all(scale))
  xscaled<-cbind(x_scaled,
        x_vars%>% dplyr::select(c(PATNO_EVENT_ID, PATNO))) 
  return(x_scaled)
}

x_vars_sel_train_lfc%>% dplyr::select(c(PATNO_EVENT_ID, PATNO))
x_vars_sel_train_scaled<-scale_numeric(x_vars_sel_train_lfc)

x_vars_sel_train_scaled_y<-cbind(x_vars_sel_train_scaled, y)
dim(x_vars_sel_train)
length(x_vars_sel_train_scaled_y$y)


x_vars_sel_test_scaled<-scale_numeric(x_vars_sel_test)
dim(x_vars_sel_test_scaled)




x_vars_sel_train_y[1:5]

#GGally::ggpairs(x_vars_sel_train_y_scaled)
# TODO: Try different models, 
# Measure error
# Try rmse
# compare to rnas 
fit <- lm(as.formula( paste('y ~ hsa.miR.101.3p + hsa.miR.126.5p+hsa.miR.144.5p+
                      hsa.miR.190a.5p+SEX+AGE+NP3_TOT+',paste(to_Add, collapse=' + ' ))),
          data = x_vars_sel_train_scaled_y)

formula_base<-'y ~SEX+AGE+NP3_TOT+con_putamen+'
fit <- lm(as.formula( paste(formula_base,paste(to_Add, collapse=' + ' ))),
          data = x_vars_sel_train_scaled_y)

print(formula_base)







svmfit = svm(y ~ ., data = dat, kernel = "linear", cost = 10, scale = FALSE)

fit <- svm(as.formula( paste(formula_base,paste(to_Add, collapse=' + ' ))),
                                  data = x_vars_sel_train_scaled_y)



formula_base<-'y ~SEX+AGE+NP3_TOT+con_putamen'
fit <- lm(as.formula( paste(formula_base)),
          data = x_vars_sel_train_scaled_y)
library(e1071)

x_vars_sel_test_scaled$PATNO

predicted<-stats::predict(fit, x_vars_sel_test_scaled)
predicted
y_test

x_vars_sel_test_scaled$
cbind(y_test, predicted)



error<-sqrt(mean((y_test - predicted)^2,na.rm=TRUE))
error


fit$coefficients
fit_an<-anova(fit)
fit_an
f <- summary(fit)
f$r.squared
f



x_vars_df<-as.data.frame(x_vars_sel_train_scaled)
x_vars_df

x_vars_df[,1:5] %>% mutate_all(as.numeric)


hist(x_vars_df[,1:5] %>% mutate_all(as.numeric))



hist(x_vars_df[,4])




