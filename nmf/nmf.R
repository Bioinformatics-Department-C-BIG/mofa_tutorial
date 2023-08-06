
setwd('/tmp/')
library(FactoMineR)
library(pixmap)
install.packages('NMF')
library(NMF)


# TODO: setup the home dir 
# TODO save results 
# TODO: make plots 
mod='RNA'
mod='miRNA'
mod='proteomics'

x1_se<-mofa_multi[, , mod]

x1=assays(x1_se)[[mod]]


mofa_multi

mofa_multi$COHORT

r.mfa <- FactoMineR::MFA(x1)



res<-NMF::nmf(x1, 8)
outdir
## MULTI RUN
NFACTORS=3
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

h<-coef(res)
h[1,]
x1_se$PATNO


### Correlations for each factor
apply(h, 1, cor, x=x1_se$COHORT)

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






### enrichment analysis? 







