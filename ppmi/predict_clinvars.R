
library(GGally)
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
selected_feats2<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.995),])
selected_feats2

x_prots<-MOFAobject@data$miRNA[[1]][selected_feats2, ]

x_prots
x_feats<-rbind(x_rnas, x_prots)
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

source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
se_filt_V06<-filter_se(se, VISIT='V06', sel_coh,sel_ps)
se_filt_V04<-filter_se(se, VISIT='V04', sel_coh,sel_ps)

### Selected features from baseline
# TODO: split by caret package 
library(caret)
set.seed(3456)
# split so that there are equal numbers of different stages of NHY 
# For categorical use cut function
common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )

se_filt_V08_pd=preprocess_visit(se_filt_V08, common=common,feat_names = selected_feats)
trainIndex <- createDataPartition(se_filt_V08_pd$NHY, p = .8, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)  
  
x_sel_bl<-se_filt_BL[selected_feats, training_samples]

x_vars<-cbind(t(MOFAobject@data$RNA[[1]][selected_feats, ]))
x_vars<-cbind(t(x_feats), MOFAobject@samples_metadata[, c('AGE')] )

x_vars<-cbind(t(x_feats), MOFAobject@samples_metadata[, c('AGE')] )




#### Outcome is from Visit 8
#### 1. V
#
#
# TODO: train test split: 
# 1. subset by samples, and by features
# ALSO split to train test
y= se_filt_V08[]

dim(x_vars)
dim(x_vars)

y
fit <- lm(y ~ Z)
fit <- lm(y ~ x_vars)

fit$coefficients
fit_an<-anova(fit)
fit_an
f <- summary(fit)
f$r.squared
f

sel_factors
