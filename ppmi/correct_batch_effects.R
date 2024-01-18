#BiocManager::install('sva')

library('sva')
library('WGCNA')
library('R.filesets')

library('FactoMineR')


# https://docplayer.net/21369417-Differential-analysis-of-count-data-the-deseq2-package.html
### Analysis of batch effects  due to ###
# 1. remove the most variable genes from the pools (filtered_genes)
# Correct for sequencing metrics: 
# 1. usable bases, 2. plate 
# (3. site is also associated witht echnical variability but also biological so we don't remove it..###)


# remove genes related to the first PCA ###

filtered_genes<-read.csv(paste0(data_dir, 'ppmi/ppmi_data/rnaseq/filteredGenes.csv'))
remove_genes<-filtered_genes$perc0.1


cell_corr=FALSE
prefix='rnas_'; process_mirnas<-FALSE;  
#prefix='mirnas_'; process_mirnas<-TRUE;  
if (cell_corr){
    to_remove_covars<-c('Usable_Bases_SCALE', 'Plate', 'Neutrophil.Score')

}else{
    to_remove_covars<-c('Usable_Bases_SCALE', 'Plate')

}
to_keep_covars = c('AGE_SCALED', 'SEX', 'COHORT')


source(paste0(script_dir,'ppmi/config.R'))

input_file_mirs_norm
se_mirs_prenorm = load_se_all_visits(input_file = input_file_mirs_norm, combined=combined_bl_log)
se_pr_mirs_prenorm<-preprocess_se_deseq2(se_mirs_prenorm, min.count = min.count) # scale and transform covariates used in preprocessing 

input_file
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log)
# TODO: here try also the tpm measures that are already normalized!! 
se_pr <- preprocess_se_deseq2(se, min.count = 20) # scale and transform covariates used in preprocessing 
se_pr <- preprocess_se_deseq2(se, min.count = min.count) # scale and transform covariates used in preprocessing 




###
se_pr$Usable_Bases_SCALE<-as.numeric(se_pr$Usable_Bases_SCALE)

ddsSE <- DESeqDataSet(se_pr, 
                      design = as.formula(formula_deseq))
ddsSE<-estimateSizeFactors(ddsSE)


######
batch<-ddsSE$Usable_Bases_SCALE
cohorts<-ddsSE$COHORT

removedCovariates<-colData(ddsSE)[,to_remove_covars]
retainedCovariates<-colData(ddsSE)[,to_keep_covars]




## PCA with 1. vsd 2. CPM 3. 
if (file.exists(vst_cor_all_vis)){
  
  # load the already computed files 

   vsd_cor<-loadRDS(vst_cor_all_vis)
   vsd<-loadRDS( vst_all_vis )

}else{

    # TODO: turn into function to create the correction 
    if (NROW(assay(ddsSE)) > 1000){
       vsd<-vst(ddsSE)# switched tot this because it is faster 

    }else{
      vsd<-varianceStabilizingTransformation(ddsSE,blind = FALSE)
    }
    saveRDS(vsd,  vst_all_vis )


    # correction 
    adjusted_data<-empiricalBayesLM(t(as.matrix(assay(vsd))), removedCovariates=removedCovariates, 
                                    retainedCovariates = retainedCovariates )
    vsd_cor<-vsd
    assay(vsd_cor)<-t(adjusted_data$adjustedData)

    ## Save the vsd corrected with all visits inside 
  

  #### if you want to skip re-running load 
    if (!process_mirnas){
      vsd_cor_filt<-vsd_cor[!(rownames(vsd_cor) %in% filtered_genes$perc0.1),]
      #
      #vsd_cor_filt<-vsd_cor_filt[rownames(vsd_cor_filt)]    

      dim(vsd_cor_filt)
      
     
      vsd_cor<-vsd_cor_filt
      
      assay_r<-gsub( '\\..*','' ,rownames(assay(vsd_cor)))
      vsd_cor<-se_filt[assay_r %in% intersect(assay_r, remaining_genes),]



    } 
    saveRDS(vsd_cor, vst_cor_all_vis)

}


### filter by visit
# 1.  
# Compare corrected and uncorected only for V08
vsd_V08<-filter_se(se = vsd,VISIT = 'V08',sel_coh = sel_coh, sel_sub_coh = sel_subcoh) 
vsd_cor_V08<-filter_se(se = vsd_cor,VISIT = 'V08',sel_coh = sel_coh, sel_sub_coh = sel_subcoh) 
vsd_cor_V08_rem<-vsd_cor_V08[!(rownames(vsd_cor_V08) %in% filtered_genes$perc0.1), ]


saveRDS(vsd_cor, vst_cor_all_vis_filt)

pca_files = paste0(data_dir, '/ppmi/plots/single/pca/')
pca_pars = paste0(cell_corr)



plot_corr = TRUE # Choose corrected or uncorrected plots 
if (plot_corr){
  vsd_p<-vsd_cor_V08
}else{
  vsd_p<-vsd_V08

}



      pca_data<-t(assay(vsd_p)) # transpose because rows must be samples
      meta_d<-colData(vsd_p)
      coh_d<-meta_d$COHORT_DEFINITION
      bases<-meta_d$Usable_Bases_SCALE
      pca.data_vsd <- PCA(pca_data,
                      scale.unit = TRUE, graph = FALSE)
    
      pc_all<-get_pca_ind(   pca.data_vsd )$coord
      cor(pc1,as.numeric(vsd_p$COHORT), method='spearman' )

      
      fviz_eig(pca.data_vsd, addlabels = TRUE, ylim = c(0, 70))

      graphics.off()
      covar<-'Usable_Bases_SCALE'
      covar<-'Plate'
      covar<-'Neutrophil.Score'
      covar<-'COHORT'
      covar<-'Neutrophil.Score'

      View(pca.data_vsd$loadings)
      pc_neutr_cor=cor(pc_all,as.numeric(colData(vsd_p)[,covar]), method='spearman' )
  
      pc_ind_p<-fviz_pca_ind(pca.data_vsd,col.ind = as.numeric(colData(vsd_p)[,covar]), 
                             label=covar
                              )+
        labs(title = paste0("PCA, color: ", covar, '\ncor: ', paste0(round(pc_neutr_cor[1:2], digits=2 ),  collapse=', ')) )      
      
      
      pc_ind_p
      ##pc_ind_ps[[i]]=pc_ind_p
      ggsave(paste0(pca_files, 'individuals_correct_',plot_corr, pca_pars, covar,'.jpeg'), width=4, height = 4, dpi=300)
      
      ### graph of variables 
      
      
      
      fviz_pca_var(pca.data_vsd,
        col.var="contrib")+
        scale_color_gradient2(low="blue", mid="white", 
                              high="red", midpoint=1.5)+theme_bw()
      
      
      #### TODO: calculate PCs for ALL the datasets 
      
      
      ggsave(paste0(pca_files, 'variables', pca_pars,'.jpeg'), width=6, height = 6)
      
      
      


# todo: load corrected without cell types and corrected with cell types 

vsd1=vsd_V08
vsd2=vsd_cor_V08

####input mirs alrready normalized
se_pr_mirs_prenorm_log<-se_pr_mirs_prenorm
assay(se_pr_mirs_prenorm_log)=log2(assay(se_pr_mirs_prenorm)+1)
se_mirs_prenorm_log_V08<-filter_se(se = se_pr_mirs_prenorm_log,VISIT = 'V08',sel_coh = sel_coh, sel_sub_coh = sel_subcoh) 

vsd1=se_mirs_prenorm_log_V08;dim(vsd1)
vsd2=vsd_cor_V08; dim(vsd2)

View(assay(vsd1))
View(assay(vsd2))

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


## Choose df, variance or mean
df<-data.frame(vsd=median_expr_vsd,vsd_corrected=median_expr_vsd_cor)
dim(df)
df<-cbind(df, colData(vsd_p)[, c('Usable_Bases_SCALE', 'Plate')] );
colnames(df)
df_melt<-reshape2::melt(df,id.vars=c('Usable_Bases_SCALE', 'Plate') )
colnames(df_melt)


which(df_melt$value <8)
df_melt_remove_outl<-df_melt[df_melt$value <10 &df_melt$value >8 ,]
ggplot(df_melt_remove_outl, aes(y=value, fill=Plate))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 2)

ggsave(pca_files, '/var.csv')

ggplot(df_melt, aes(y=value, fill=Plate))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 2)

ggplot(df_melt_remove_outl, aes(y=value, fill=Plate))+
  geom_boxplot()


ggplot(df_melt, aes(y=value, x=Usable_Bases_SCALE))+
  geom_point()+
  
  facet_wrap(~variable, nrow = 2)




########### CHECK POOLS per plate ###########
# 1. Do the pools on each plate have different distributions of gene expression 
# Could not find pool samples 1009 and 1010--> 
# TODO: 1. check the raw data for 1009 and 1010

colData(se[,se$PATNO %in% c(1009, 1010)])




### METADATA correlations 
## 1. 
clinical_data$Usable_Bases_SCALE
selected_covars2_progression[!selected_covars2_progression %in% colnames(colData(se))]


se$Usable_Bases_SCALE<-as.numeric(scale(se$`Usable.Bases....`))
selected_covars2_cor<-c('SITE', 'Plate', 'NP2PTOT', 'MCATOT', 'Neutrophil.Score', 'Usable_Bases_SCALE')
clinical_data<-as.data.frame(colData(se)[, selected_covars2_cor])

clinical_data$Plate<-as.factor(clinical_data$Plate)
#clinical_data = as.data.frame(sapply((clinical_data), as.factor))
clinical_data = sapply((clinical_data), as.numeric)
clinical_data
colnames(clinical_data)



cor <- psych::corr.test(clinical_data, clinical_data, method = "pearson", adjust = "BH" )
stat <- cor$r

jpeg(paste0(outdir, '/covariates/corr_plot.jpeg'))
cor_plot<-corrplot::corrplot(stat, tl.col = "black", title="Pearson correlation coefficient",diag=FALSE, )

dev.off()


hist(se$Neutrophil.Score)











































































