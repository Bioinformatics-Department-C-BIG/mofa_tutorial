
# https://bioconductor.org/packages/devel/bioc/vignettes/granulator/inst/doc/granulator.html
# estimation of cell types with granulator using rna seq 

#BiocManager::install('granulator')
library(granulator)
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
  ABIS_S1 = sigMatrix_ABIS_S1, 
  ABIS_S2 = sigMatrix_ABIS_S2, 
  ABIS_S3 = sigMatrix_ABIS_S3)



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


# deconvolute input data using all available methods by default
decon <- deconvolute(m = bulkRNAseq, sigMatrix = sigList, methods=get_decon_methods()[1:4])


qprog_ABIS_S0



library('psych')

### apply to the estimated 
names(decon$proportions)
for (estimated in decon$proportions ){

      #estimated<-decon$proportions$nnls_ABIS_S0
      colnames(estimated)
      max(decon$proportions$nnls_ABIS_S0)

      # plot cell type proportions for svr model on ABIS_S0 reference profile
      plot_proportions(deconvoluted = decon, method = 'nnls', signature = 'ABIS_S3')
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






x_covars$age
x_covars$age_at_visit
fit_lm<-lm(y~x_covars$age_at_visit+x_covars$SEX+x_covars$Usable.Bases....+x+x_covars$EVENT_ID+x)


y
#fit_lm<-glm(y~x_covars$age_at_visit+x_covars$SEX+x_covars$Usable.Bases....+x+x_covars$EVENT_ID+x,  family = "binomial")


fit_lm

plot(y ~ x)

#correlate_factors_with_covariates
sm[variates_to_p]
variates_to_p_cors<-c(variates_to_p, 'Multimapped....', 'Uniquely.mapped....', 'SITE', 'RIN.value', 
    'Usable_bases_SCALE')
cor <- psych::corr.test(sm[variates_to_p], sm[variates_to_p], method = "pearson", 
        adjust = "BH")
stat<-cor$r
stat
 corrplot(stat, tl.col = "black", title = "Pearson correlation coefficient")















library(corrplot)
M=cor(sm[variates_to_p])
corrplot(sm[variates_to_p])










