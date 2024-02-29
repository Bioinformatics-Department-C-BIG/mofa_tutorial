
#install.packages('psych')
#####  MOFA ANALYSIS 

### this one depends on mofa application to inherit: 
### 1. MOFAobject, 2. factors, 
# 3. clinical variables
# 4. outdirs 
# 5. csf/plasma/untargeted flags
#source('enrichment.R')


#### MOFA ANALYSIS ####
### 1. get correlations with covariates 

###### CORRELATIONS WITH COVARIATES ##########
##### Corelations
## TODO: load mofa object from outdir 

library('pheatmap')
#library('kml')
library('dplyr')
library('factoextra')


#### Covariance of factors with metadata 
#source(paste0(script_dir,'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir,'ppmi/mofa_utils.R'))
source(paste0(script_dir,'ppmi/plotting_utils.R'))

outdir

vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[1]]

write_vars_output(MOFAobject, vars_by_factor) # output plots on var metrics


### ADD FUTURE CHANGES ####
############################
sm<-samples_metadata(MOFAobject)
#res<-df_change[match(sm$PATNO, df_change$PATNO),]
patno_event_ids=samples_metadata(MOFAobject)$PATNO_EVENT_ID

sm<-samples_metadata(MOFAobject);

# diff variables holds the vars with the future changes 


#### TODO: add clinical clusters as a function 
#1. 

#MOFAobjectBL <- subset_groups(MOFAobject, groups = 1)
# PREPROCESS

MOFAobject@samples_metadata$Biological.group<-factor(MOFAobject@samples_metadata$Biological.group)

sm$Biological.group
patient_groups<-c('TF', 'TnF')
Patient_samples_only<-MOFAobject@samples_metadata$sample[MOFAobject@samples_metadata$Biological.group %in%patient_groups ]
Patient_samples_only



MOFAobjectPatients <- subset_samples(MOFAobject, samples=Patient_samples_only)



sm<-samples_metadata(MOFAobject)
sm






correlate_factors_categorical<-function(y, MOFAobject=MOFAobject){
  ##  obtain correlations of factors with the clinical clusters 
  #' explicitly done for categorical variables
  #' @param name 
  #'
  #'
  #'
  #'
    clust_name=paste0(y)
    #clust_name='Biological.group'
 
      sm<-samples_metadata(MOFAobject)
      if (clust_name %in% colnames(sm)){
        covariates=samples_metadata(MOFAobject)[,clust_name ]
        Z <- get_factors(MOFAobject, factors = 1:5, 
                         as.data.frame = FALSE)

        Z <- do.call(rbind, Z)

        #psych::corr.test(Z, covariates, method = "pearson", 
        #                  adjust = "BH")
        y_pvals<-apply(Z, 2, function(z) {
          kruskal.test(z,covariates)$p.value })
        y_pvals<-as.numeric(p.adjust(y_pvals, method = 'BH'))
        #y_pvals<-y_pvals<0.05
        # do Kruskal Wallis test to see whether or not there is statistically significant difference between three or more groups

        return(y_pvals)
      }
  }
  

# TODO: padjust
clinical_scales=c('NP2PTOT', 'NP3TOT', 'SCAU_TOT', 'RBD_TOT', 'updrs3_score', 'updrs2_score')
add_clinical_clusters=FALSE

if (add_clinical_clusters){
  
  cors_kruskal<-sapply(clinical_scales, correlate_factors_categorical)
  rownames(cors_kruskal)<-1:N_FACTORS
  kw_cors<-format(cors_kruskal, digits=1)<0.05
  factors_cor_with_clusters<-kw_cors
  #library(pheatmap)
  pheatmap(-log10(cors_kruskal))
  #}
  
  #)
}






# CORS PATIENTS ONLY #### 
dir.create(paste0(outdir, '/covariates/'))
MOFAobjectPD_sel = MOFAobject
MOFAobject_sel=MOFAobject
stats<-apply(MOFAobjectPD_sel@samples_metadata, 2,table )
non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )
covariates_dir
covariates_dir<-paste0(outdir, '/covariates/covariate_corelations_')
covars_f_pearson<-paste0(covariates_dir, 'pearson.csv' )
covars_f_pvalue<-paste0(covariates_dir, 'pvalue.csv' )
covars_f_pearson_pd<-paste0(covariates_dir, 'pearson_pd.csv' )
covars_f_pvalue_pd<-paste0(covariates_dir, 'pvalue_pd.csv')

 names(non_na_vars)
 force_cors=TRUE
####TODO: maybe filter out some clinvars or take the most important because it takes a while....!#####
if (file.exists(covars_f_pearson_pd) & !(force_cors)){
  # Loading covariates from file
  print('Load covariates from file')
  cors_pearson_pd<-read.csv2(covars_f_pearson_pd, row.names=1)
  cors_all_pd<-read.csv2(covars_f_pvalue_pd, row.names=1)
  cors<-read.csv2(covars_f_pvalue, row.names=1)
  cors_pearson<-read.csv2(covars_f_pearson, row.names=1)
  cors_pearson
}else{
  stats<-apply(MOFAobjectPD_sel@samples_metadata, 2,table )
  sm_pd<-MOFAobjectPD_sel@samples_metadata
  # Checks if the variance >0 to obtain corelation
  non_na_vars<-which(!is.na(sapply(stats,mean)) & apply(sm_pd,2,var, na.rm=TRUE)>0 ) 
   length(names(non_na_vars) )

  cors_both<-get_correlations(MOFAobjectPD_sel, names(non_na_vars))
  cors_pearson_pd = as.data.frame(cors_both[[2]]); # H
  cors_all_pd = as.data.frame(cors_both[[1]])


  write.csv2(cors_pearson_pd, covars_f_pearson_pd)
  write.csv2(cors_all_pd, covars_f_pvalue_pd )

# CORS ALL SAMPLES 
  stats<-apply(MOFAobject_sel@samples_metadata, 2,table )
  sm<-MOFAobject_sel@samples_metadata
  non_na_vars<-which(!is.na(sapply(stats,mean)) & apply(sm,2,var, na.rm=TRUE)>0 )
  cors_both<-get_correlations(MOFAobject_sel, names(non_na_vars))
  cors_both[[2]]
  cors_pearson=as.data.frame(cors_both[[2]]); cors=as.data.frame(cors_both[[1]]); cors_all=cors_both[[1]]

  write.csv2(cors_pearson,covars_f_pearson)
  write.csv2(cors,covars_f_pvalue)

}

################ DEFINE SETS OF FACTORS WITH DIFFERENT NAMES  ####
# corelate categorical covariates


















