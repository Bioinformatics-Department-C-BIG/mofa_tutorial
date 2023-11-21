#BiocManager::install('sva')

library('sva')
library('WGCNA')
library('R.filesets')

library('FactoMineR')
### Analysis of batch effects  due to ###
# 1. remove the most variable genes from the pools (filtered_genes)
# Correct for sequencing metrics: 
# 1. usable bases, 2. plate 
# (3. site is also associated witht echnical variability but also biological so we don't remove it..###)
# remove genes related to the first PCA ###

filtered_genes<-read.csv(paste0(data_dir, 'ppmi/ppmi_data/rnaseq/filteredGenes.csv'))
remove_genes<-filtered_genes$perc0.1

formula_deseq<-"~AGE_SCALED+SEX+Usable_Bases_SCALE+COHORT"



prefix='mirnas'
input_file=mirnas_all_visits_fname
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log)

# TODO: here try also the tpm measures that are already normalized!! 
se_pr<-preprocess_se_deseq2(se, min.count = 20) # scale and transform covariates used in preprocessing 
se_pr<-preprocess_se_deseq2(se, min.count = MIN_COUNT_M) # scale and transform covariates used in preprocessing 



#kruskal.test(se_pr$Usable_Bases_SCALE, as.factor(se_pr$Plate))

se_pr$Usable_Bases_SCALE<-as.numeric(se_pr$Usable_Bases_SCALE)

ddsSE <- DESeqDataSet(se_pr, 
                      design = as.formula(formula_deseq))
ddsSE<-estimateSizeFactors(ddsSE)


######
batch<-ddsSE$Usable_Bases_SCALE
cohorts<-ddsSE$COHORT
retainedCovariates<-colData(ddsSE)[,c('AGE_SCALED', 'SEX', 'COHORT')]
removedCovariates<-colData(ddsSE)[,c('Usable_Bases_SCALE', 'Plate')]

## PCA with 1. vsd 2. CPM 3. 

vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
### Asjustment works on gaussian data so insert vsd or log cplm
adjusted_data<-empiricalBayesLM(t(as.matrix(assay(vsd))), removedCovariates=removedCovariates, 
                                retainedCovariates = retainedCovariates )

vsd_cor<-vsd
assay(vsd_cor)<-t(adjusted_data$adjustedData)

## Save the vsd corrected with all visits inside 

saveRDS(vsd_cor,vst_cor_all_vis)
vsd_cor<-loadRDS(vst_cor_all_vis)
dim(vsd_cor)
vsd_cor_filt<-vsd_cor[!(rownames(vsd_cor) %in% filtered_genes$perc0.1),]
saveRDS(vsd_cor_filt,vst_cor_all_vis_filt)

### filter by visit
# 1. 
# Compare corrected and uncorected only for V08
vsd_V08<-filter_se(se = vsd,VISIT = 'V08',sel_coh = sel_coh, sel_sub_coh = sel_subcoh) 
vsd_cor_V08<-filter_se(se = vsd_cor,VISIT = 'V08',sel_coh = sel_coh, sel_sub_coh = sel_subcoh) 
vsd_cor_V08_rem<-vsd_cor_V08[!(rownames(vsd_cor_V08) %in% filtered_genes$perc0.1), ]
saveRDS(vsd_cor,vst_cor_all_vis_filt)



vsd_p<-vsd_cor_V08

      pca_data<-t(assay(vsd_p)) # transpose because rows must be samples
      meta_d<-colData(vsd_p)
      coh_d<-meta_d$COHORT_DEFINITION
      bases<-meta_d$Usable_Bases_SCALE
      pca.data_vsd <- PCA(pca_data,
                      scale.unit = TRUE, graph = FALSE)
    
     
      
      fviz_eig(pca.data_vsd, addlabels = TRUE, ylim = c(0, 70))

      graphics.off()
      covar<-'Usable_Bases_SCALE'
      covar<-'Plate'
      
      pc_ind_p<-fviz_pca_ind(pca.data_vsd,col.ind = as.numeric(colData(vsd_p)[,covar]), 
                             label=covar
                              )+
        labs(title = paste0("PCA, color: ", covar))
      
      
      
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
      
      
      

vsd1=vsd_cor_V08
vsd2=vsd_cor_V08_rem
median_expr_vsd_cor<-apply(t(assay(vsd1)), 1, mean)
median_expr_vsd<-apply(t(assay(vsd2)), 1, mean)

sd_vsd<-apply(t(assay(vsd1)), 1, sd)
sd_vsd_cor<-apply(t(assay(vsd2)), 1, sd)
length(sd_vsd)
length(sd_vsd_cor)

### Plot the difference in mean and sd between the raw and corrected datasets 
##
df<-data.frame(vsd=sd_vsd, vsd_corrected=sd_vsd_cor)
df<-data.frame(vsd=median_expr_vsd, vsd_corrected=median_expr_vsd_cor)
df<-data.frame(vsd=median_expr_vsd,vsd_corrected=median_expr_vsd_cor)
dim(df)
df<-cbind(df, colData(vsd_p)[, c('Usable_Bases_SCALE', 'Plate')] );
colnames(df)
df_melt<-reshape2::melt(df,id.vars=c('Usable_Bases_SCALE', 'Plate') )
colnames(df_melt)

ggplot(df_melt, aes(y=value, fill=Plate))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 2)



ggplot(df_melt, aes(y=value, x=Usable_Bases_SCALE))+
  geom_point()+
  
  facet_wrap(~variable, nrow = 2)




########### CHECK POOLS per plate ###########
# 1. Do the pools on each plate have different distributions of gene expression 
# Could not find pool samples 1009 and 1010--> 
# TODO: 1. check the raw data for 1009 and 1010

colData(se[,se$PATNO %in% c(1009, 1010)])


