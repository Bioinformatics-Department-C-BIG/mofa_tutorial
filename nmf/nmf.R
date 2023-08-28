
setwd('/tmp/')
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
mod='RNA'


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
NFACTORS=8
nrun=5


## SAVE AND LOAD 
out_nmf<-paste0(outdir,'/../','multirun_', NFACTORS, '_', nrun,'_', 
                    mod)




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



colnames(cor$r) 
rownames(cor$sef)
cors_non_na<-names(which(!is.na( cor$p[,1])))
T<--log10(0.00)
sig<-which( rowMins(cor$p)<0.05 & rowMaxs(cor$r)<0.95 ) 
which( rowMins(cor$p)<0.05)

cor$r['CONCOHORT',]
stat[sig,]


plot='r'


nmf_param_str<-paste0('nmf/plots/','cor_', mod, '_', NFACTORS )
if (plot=="r") {
  stat <- cor$r
  png(paste0( nmf_param_str, '.png' ))
  
  corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")
  dev.off()  
  
  write.csv(stat[sig,], paste0('nmf/plots/cor_',mod, '.csv'  ))
} else if (plot=="log_pval") {
  stat <- cor$p
  stat[stat>alpha] <- 1.0
  if (all(stat==1.0)) stop("All p-values are 1.0, cannot plot the histogram")
  stat <- -log10(stat)
  stat[is.infinite(stat)] <- 1000
  if (transpose) stat <- t(stat)
  if (return_data) return(stat)
  col <- colorRampPalette(c("lightgrey", "red"))(n=100)
  pheatmap::pheatmap(stat, main="log10 adjusted p-values", cluster_rows = FALSE, color=col, ...)
  
} else {
  stop("'plot' argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
}




## Get back a table for all metadata ? 


### data matrix 

V.hat<-fitted(res)
dim(V.hat)

summary(res)


###prior knowledge on class labels
summary(res,class=x1_se$COHORT)


w<-basis(res)
w
### 


s<-featureScore(res)

s




### enrichment analysis? 







