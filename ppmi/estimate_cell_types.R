
# https://bioconductor.org/packages/devel/bioc/vignettes/granulator/inst/doc/granulator.html
# estimation of cell types with granulator using rna seq 

#BiocManager::install('granulator')

library(granulator)
library(R.filesets)
script_dir<-paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../../')
script_dir
cell_corr_deseq = TRUE
#VISIT='V08'
source(paste0('ppmi/setup_os.R'))
#source(paste0(script_dir, 'ppmi/setup_os.R'))
script_dir
source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
# load datasets for deconvolution of PBMC RNA-seq data
load_ABIS()


# print TPM values in bulk RNA-seq
bulkRNAseq_ABIS[1:5, 1:5]

# print TPM values in reference profile matrix
sigMatrix_ABIS_S0[1:5, 1:5]


# print measured cell type proportions (percentages)
groundTruth_ABIS[1:5, 1:5]



# create list if multiple signature matrices to test simultaneously
sigList = list(
  ABIS_S0 = sigMatrix_ABIS_S0,
  ABIS_S1 = sigMatrix_ABIS_S1
  #ABIS_S2 = sigMatrix_ABIS_S2, 
  #ABIS_S3 = sigMatrix_ABIS_S3)
)


  # plot signature matrix similarity matrices
plot_similarity(sigMatrix=sigList)

vst_cor_all_vis_data<-loadRDS(vst_cor_all_vis)


process_mirnas=FALSE; 
source(paste0(script_dir, '/ppmi/config.R'))

se_rnas<-load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_rnas_cpm<-cpm(assay(se_rnas))

se_rnas_cpm
##
bulkRNAseq<-assay(se_rnas)
bulkRNAseq<-se_rnas_cpm

rownames(bulkRNAseq)
rownames(bulkRNAseq)<-get_symbols_vector(rownames(bulkRNAseq))

get_decon_methods()[4:7]

# deconvolute input data using all available methods by default
decon_file<-paste0(output_files, '/decon.RDS')
if (file.exists(decon_file)){
  decon<-loadRDS(decon_file)}else{
  decon <- deconvolute(m = bulkRNAseq, sigMatrix = sigList, methods=get_decon_methods()[4:7])

  saveRDS(decon, decon_file)
}




library('psych')

### apply to the estimated 
names(decon$proportions)
for (estimated in decon$proportions ){

      #estimated<-decon$proportions$nnls_ABIS_S0
      colnames(estimated)
      max(decon$proportions$nnls_ABIS_S0)

      # plot cell type proportions for svr model on ABIS_S0 reference profile
    #  plot_proportions(deconvoluted = decon, method = 'nnls', signature = 'ABIS_S3')
      # plot cell type proportions
      plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)

      ground_truth<-as.data.frame(colData(se_rnas))

      ground_truth$PATNO_EVENT_ID
      match(rownames(estimated),ground_truth$PATNO_EVENT_ID )


      match(ground_truth$PATNO_EVENT_ID , rownames(estimated))


      estim_matched<-estimated[match(ground_truth$PATNO_EVENT_ID , rownames(estimated)),]

      estim_matched$Neutrophils.LD

      se_rnas$Neutrophil.Score

      cor_neu<-corr.test(estim_matched$Neutrophils.LD, ground_truth$`Neutrophils....`)
      
      corr.test(estim_matched$Neutrophils.LD, ground_truth$`Neutrophils....`)$p
      print(paste(cor_neu$r, cor_neu$p))

}

sm<-samples_metadata(MOFAobject)
sm<-samples_metadata(MOFAobject_clusts)



estimated<-decon$proportions$qprogwc_ABIS_S0

estim_matched<-estimated[match(ground_truth$PATNO_EVENT_ID , rownames(estimated)),]

corr.test(estim_matched$Neutrophils.LD, ground_truth$Neutrophil.Score, )

corr.test(estim_matched$Neutrophils.LD, ground_truth$Neutrophil.Score)$r
corr.test(estim_matched$Neutrophils.LD, ground_truth$Neutrophil.Score)$p.adj


ground_truth$INEXPAGE
filt_samples<-ground_truth$COHORT%in% c(1,2) & ground_truth$EVENT_ID%in% c('BL','V04', 'V06','V08') & (ground_truth$INEXPAGE %in% c('INEXPD', 'INEXHC') )


filt_samples
cestim_matched
cor_cohort <- corr.test(estim_matched[filt_samples,], ground_truth$COHORT[filt_samples])
cor_cohort$p.adj<0.05
cor_cohort$r[cor_cohort$p.adj<0.05]




## regression 
y= ground_truth$COHORT[filt_samples]
x_covars<-ground_truth[filt_samples, ]

i=8
colnames(estim_matched)
colnames(estim_matched)[i]
x=estim_matched[filt_samples,i]
x

fit_lm<-lm(y~x_covars$age_at_visit+x_covars$SEX+x_covars$Usable.Bases....+x+x_covars$EVENT_ID+x)


#correlate_factors_with_covariates
# List median cell type estimation by disease control 

# Correlate 
# Estimate and check for correlations or significance between the two groups 
covariates_to_p<-c('AGE', 'SEX', 'NP2PTOT', 'NHY', '')

variates_to_p<-colnames(estimations)
variates_to_p_cors<-c(variates_to_p, 'Multimapped....', 'Uniquely.mapped....',  'RIN.Value', 
    'Usable_Bases_SCALE')

   variates_to_p_cors %in%  colnames(sm)
 sm[variates_to_p_cors]==0
cor <- psych::corr.test(sm[variates_to_p_cors], sm[variates_to_p_cors], method = "pearson", 
        adjust = "BH")
stat<-cor$r
stat
 corrplot(stat, tl.col = "black", title = "Pearson correlation coefficient")




library(corrplot)
M=cor(sm[variates_to_p])
corrplot(sm[variates_to_p])
sm2<-sm[!is.na(sm$COHORT),]
# 1. apply shapiro
# 2. Check normality 
wil1<-wilcox.test(sm$Neutrophil.Score, sm$COHORT)
wil1$p.value


wilcox_cell_types<-apply(sm2[, c(colnames(estimations), 'Neutrophils....', 'Lymphocytes....')],2 ,
function (x) {wilcox.test(x, sm2$COHORT)$p.value} )
wilcox_cell_types


variances_cells<-sm2[, c(colnames(estimations))] %>%
    summarize_all(var, na.rm=TRUE) %>%
    as.data.frame() %>% t() %>%
    as.data.frame()

    variances_cells
write.csv(round(variances_cells, digits=2), paste0(outdir, '/cell_types/', 'variances.csv'))


medians<-sm2[, c(colnames(estimations), 'COHORT', 'Neutrophils....', 'Neutrophil.Score')] %>%
    group_by(COHORT) %>%
    summarize_all(median, na.rm=TRUE) %>% t() %>%
    as.data.frame() 

    medians
write.csv(round(medians, digits=3), paste0(outdir, '/cell_types/', 'medians.csv'))
IQRS<-sm2[, c(colnames(estimations), 'COHORT', 'Neutrophils....', 'Neutrophil.Score')] %>%
    group_by(COHORT) %>%
    summarize_all(IQR, na.rm=TRUE) %>% t() %>%
    as.data.frame() 

sm2$cky
sm2<-samples_metadata(MOFAobject_sel)
sm2$NP2PTOT_LOG_clust
cell_type_df<-sm2[, c(colnames(estimations), 'COHORT', 'Neutrophils....', 'Neutrophil.Score', 'COHORT')]
cell_type_df<-sm2[, c(colnames(estimations), 'Neutrophils....', 'Neutrophil.Score', 'NP2PTOT_LOG_clust')]
cell_type_df<-sm2[, c(colnames(estimations),  'Neutrophils....', 'Neutrophil.Score','COHORT')]

y_clust_cells = 'COHORT'
cell_type_df<-cell_type_df[,!(colnames(cell_type_df) %in% c('Plasmablasts', 'T.gd.Vd2' ))]
cell_type_df_long<-melt(cell_type_df, id.vars = c(y_clust_cells))
colnames(cell_type_df_long)


cell_type_df_long

ggplot(cell_type_df_long, aes_string(group=y_clust_cells,, x='value') )+
geom_density(aes(x = value, fill=NP2PTOT_LOG_clust), alpha=0.6)+
 scale_fill_viridis_d(option='magma')+
facet_wrap(.~variable, scales='free' )
ggsave( paste0(outdir, '/cell_types/', 'hist_clusters.jpeg'))

cell_type_df_long<-cell_type_df_long[!is.na(cell_type_df_long$COHORT),]




cell_type_df_long$COHORT<-as.factor(cell_type_df_long$COHORT)



ggplot(cell_type_df_long, aes_string( y='value' ,x=y_clust_cells))+
geom_violin(aes_string(x = y_clust_cells,y='value', fill=y_clust_cells), alpha=0.9)+
 scale_fill_viridis_d(option='magma')+
 geom_pwc(aes_string(x=y_clust_cells, y='value'),method='wilcox_test', label = "p.adj.signif", 
 label.size = 2 )+
facet_wrap(.~variable, scales='free', nrow = 2 )

ggsave( paste0(outdir, '/cell_types/', 'violin_clusters', y_clust_cells, '.jpeg'), 
width=10, height=6)

# check shapiro test 



sm2$MAIT[sm2$COHORT==2]
dir.create(paste0(outdir, '/cell_types/'))
write.csv(wilcox_cell_types, paste0(outdir, '/cell_types/', 'wilcox_pd_hc.csv'))


























