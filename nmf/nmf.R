
#setwd('/tmp/')
library(FactoMineR)
library(pixmap)
#install.packages('NMF')
library(NMF)


# TODO: setup the home dir 
# TODO save results 
# TODO: make plots 
mod='miRNA'
mod='proteomics'
mod='RNA'
mod='miRNA'
mod='RNA' ; nmf_params<-g_params


#### Load the dataset ####
## get list of three mats 
data_full=prepare_multi_data(p_params, param_str_g, param_str_m, mofa_params)
# create multi experiment 
mofa_multi<-create_multi_experiment(data_full, combined_bl)

mofa_multi$COHORT


x1_se<-mofa_multi[, , mod]

x1=assays(x1_se)[[mod]]



#### select the highly variable #### 




dim(x1)
res<-NMF::nmf(x1, 8)
outdir
## MULTI RUN
NFACTORS=5
nrun=2


## SAVE AND LOAD 
out_nmf<-paste0(outdir,'/../','multirun_', NFACTORS, '_', nrun,'_', 
                    mod )

out_nmf<-paste0(outdir,'/../','multirun_', NFACTORS, '_', nrun,'_', 
                mod )


out_nmf_params<- paste0( 'g_', g_params, nmf_params, '_coh_', sel_coh_s,'_', VISIT_S)
outdir_nmf = paste0(outdir_orig, '/nmf/',out_nmf_params, '/');
dir.create(outdir_nmf)


if (file.exists(out_nmf)){
  res=loadRDS(out_nmf)
  
  
}else{
  
  res.multirun<-NMF::nmf(x1,NFACTORS,nrun=nrun )
  res=res.multirun
  saveRDS(res.multirun,out_nmf)
}
   




### return fitted model 
fit(res)

h <-as.data.frame(coef(res)) # factor coeficients for each sample 
dim(h)
h[1,]
x1_se$PATNO




#### Correlations for each factor
# 1. Tune to maximize cohort correlations ####
round(apply(h, 1, cor, x=x1_se$COHORT), 2)





covariates <- as.data.frame(lapply(colData(x1_se), as.numeric))
rownames(covariates)<-colData(x1_se)$PATNO_EVENT_ID
dim(colData(x1_se))
dim(t(covariates))

rownames(covariates)
rownames(t(h))
h_t<-as.data.frame(t(h))

library(psych)
duplicated(row.names(h_t))
duplicated(row.names(covariates))
colnames(covariates)
covariates$PATNO
dim(t(covariates))

dim(covariates)
cor <- psych::corr.test(covariates,h_t, method = "pearson", adjust = "BH")





nmf_param_str<-paste0('nmf/plots/','cor_', mod, '_', NFACTORS )
nmf_outdir<-nmf_param_str
dir.create(nmf_outdir)
dir.create(paste0(nmf_outdir, '/top_weights/'))
dir.create(paste0(nmf_outdir, '/enrichment/'))
dir.create(paste0(nmf_outdir, '/heatmaps/'))


## Get back a table for all metadata ? 


### data matrix 

V.hat<-fitted(res)
dim(V.hat)

summary(res)

#### Evaluation ####
###prior knowledge on class labels
summary(res,class=x1_se$COHORT)



### 


s<-featureScore(res)


## TUNING #### 
# Print the cophonetic coefficient , the RSS curve to decide number of factors #

estim.r <- nmf(x1, 2:6, nrun=3, seed=123456)
plot(estim.r)


## plot also annotation datasets 
#anndf<-as.data.frame(colData(x1_se)[, c('COHORT')]); rownames(anndf)=x1_se$PATNO

consensusmap(estim.r, annCol=x1_se$COHORT, labCol=NA, labRow=NA)


s
plot(estim.r)




#### Top Weights ####



### enrichment analysis? 







