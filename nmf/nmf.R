
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
g_params



dim(x1)
res<-NMF::nmf(x1, 8)
outdir
## MULTI RUN
NFACTORS=9
nrun=3


## SAVE AND LOAD 


out_nmf_params<- paste0( 'g_', g_params, nmf_params, '_coh_', sel_coh_s,'_', VISIT_S)
outdir_nmf = paste0(outdir_orig, '/nmf/',out_nmf_params, '_', NFACTORS );
out_nmf=paste0(outdir_nmf, 'model')

dir.create(outdir_nmf)


if (file.exists(outdir_nmf)){
  res=loadRDS(outdir_nmf)
  
  
}else{
  
  res.multirun<-NMF::nmf(x1,NFACTORS,nrun=nrun )
  res=res.multirun
  saveRDS(res.multirun,outdir_nmf)
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
h_t<-as.data.frame(t(h))









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
run_tuning=FALSE
if (run_tuning){
  estim.r <- nmf(x1, 2:6, nrun=2, seed=123)
  saveRDS(estim.r, './estimr')
  
  plot(estim.r)
  
  
  ## plot also annotation datasets 
  #anndf<-as.data.frame(colData(x1_se)[, c('COHORT')]); rownames(anndf)=x1_se$PATNO
  
  jpeg(paste0(out_nmf_params,'_consensus_map.jpeg'), res=300,units = 'in', width=10, height=9)
  p<-consensusmap(estim.r, annCol=x1_se$COHORT)
  dev.off()
  
  s
  plot(estim.r)
  summary(estim.r$fit,class=)
}





#### Top Weights ####



### enrichment analysis? 







