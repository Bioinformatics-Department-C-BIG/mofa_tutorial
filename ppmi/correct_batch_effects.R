#BiocManager::install('sva')

library('sva')
library('WGCNA')


### Analysis of batch effects  due to ###
# 1. usable bases, 2. plate 3. site ###
# remove genes related to the first PCA ###

formula_deseq<-"~AGE_SCALED+SEX+Usable_Bases_SCALE+COHORT"

kruskal.test(se_filt_V08_pr$Usable_Bases_SCALE, as.factor(se_filt_V08_pr$Plate))

raw_mat<-assay(se_filt)
se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
se_filt_V08_pr<-preprocess_se_deseq2(se_filt_V08)
se_filt_V08_pr$Usable_Bases_SCALE<-as.numeric(se_filt_V08_pr$Usable_Bases_SCALE)

ddsSE <- DESeqDataSet(se_filt_V08_pr, 
                      design = as.formula(formula_deseq))
ddsSE<-estimateSizeFactors(ddsSE)


######
batch<-ddsSE$Usable_Bases_SCALE
cohorts<-ddsSE$COHORT
retainedCovariates<-colData(ddsSE)[,c('COHORT')]
removedCovariates<-colData(ddsSE)[,c('Usable_Bases_SCALE', 'Plate')]



ddsSE_cpm<-cpm(assay(ddsSE))


## PCA with 1. vsd 2. CPM 3. 

vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
### Asjustment works on gaussian data so insert vsd or log cplm
adjusted_data<-empiricalBayesLM(t(as.matrix(assay(vsd))), removedCovariates=removedCovariates, 
                                retainedCovariates = retainedCovariates )

vsd_cor<-adjusted_data$adjustedData
dim(vsd_cor)
dim(vsd)
      pca_data<-t(assay(vsd)) # transpose because rows must be samples
      pca_data<-vsd_cor # transpose because rows must be samples
      
      meta_d<-colData(se_filt_V08_pr)
      coh_d<-meta_d$COHORT_DEFINITION
      bases<-meta_d$Usable_Bases_SCALE
      pca.data_vsd <- PCA(pca_data,
                      scale.unit = TRUE, graph = FALSE)
      
      pca.data_vsd <- PCA(pca_data,
                          scale.unit = TRUE, graph = FALSE)
     
      
      fviz_eig(pca.data_vsd, addlabels = TRUE, ylim = c(0, 70))

      graphics.off()
      covar<-'Usable_Bases_SCALE'
      pc_ind_p<-fviz_pca_ind(pca.data_vsd,col.ind = as.numeric(colData(vsd)[,covar]), 
                             label=covar
                              )
      
      
      pc_ind_p
      pc_ind_ps[[i]]=pc_ind_p
      ggsave(paste0(pca_files, 'individuals', pca_pars, '.jpeg'), width=6, height = 6)
      
      ### graph of variables 
      
      
      
      fviz_pca_var(pca.data,
        col.var="contrib")+
        scale_color_gradient2(low="blue", mid="white", 
                              high="red", midpoint=1.5)+theme_bw()
      
      
      #### TODO: calculate PCs for ALL the datasets 
      
      
      ggsave(paste0(pca_files, 'variables', pca_pars,'.jpeg'), width=6, height = 6)
}


median_expr_vsd_cor<-apply(vsd_cor, 1, mean)
median_expr_vsd<-apply(t(assay(vsd)), 1, mean)
df<-data.frame(vsd=median_expr_vsd, vsd_corrected=median_expr_vsd_cor)
df<-cbind(df, removedCovariates );
colnames(df)
df_melt<-reshape2::melt(df,id.vars=c('Usable_Bases_SCALE', 'Plate') )
colnames(df_melt)

ggplot(df_melt, aes(y=value, fill=Plate))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 2)



ggplot(df_melt, aes(y=value, x=Usable_Bases_SCALE))+
  geom_point()+
  
  facet_wrap(~variable, nrow = 2)






