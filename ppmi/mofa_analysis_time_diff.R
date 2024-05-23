
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

clinical_scales<-scale_vars_diff

MOFAobject@samples_metadata$Usable_Bases_SCALE<-scale(MOFAobject@samples_metadata$`Usable.Bases....`)
samples_metadata(MOFAobject)$SCAU26CT<-as.factor(tolower(samples_metadata(MOFAobject)$SCAU26CT))
#samples_metadata(MOFAobject)$months<-unlist(EVENT_MAP[samples_metadata(MOFAobject)$EVENT_ID], use.names = FALSE)

write.csv(round(vars_by_factor_all[[1]][[1]], digits=2), paste0(outdir,'total_variance.csv'))

vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[1]]


plot_data_overview(MOFAobject)
ggsave(paste0(outdir,'/data_overview.jpeg'))
write_vars_output(MOFAobject, vars_by_factor) # output plots on var metrics


### ADD FUTURE CHANGES ####
############################
sm<-samples_metadata(MOFAobject)
#res<-df_change[match(sm$PATNO, df_change$PATNO),]
patno_event_ids=samples_metadata(MOFAobject)$PATNO_EVENT_ID
df_change_total<-get_diff_zscores(patno_event_ids = samples_metadata(MOFAobject)$PATNO_EVENT_ID,imaging_variables_diff =imaging_variables_diff, 
                                  scale_vars_diff=scale_vars_diff )
samples_metadata(MOFAobject) = cbind(sm,df_change_total )
sm<-samples_metadata(MOFAobject);

# diff variables holds the vars with the future changes 
diff_variables<-colnames(sm)[grep('diff', colnames(sm) )]


#### TODO: add clinical clusters as a function 
#1. 
sm<-samples_metadata(MOFAobject)
y=clinical_scales[10]
y='NP2PTOT'

#sapply(clinical_scales, function(y){
#  y
# TODO: loop
all_event_ids_p<-c('BL','V02','V04','V06','V08','V10','V12','V14','V16', 'V18')
sm=MOFAobject@samples_metadata
### obtain all patient event ids to getr one row per patient!! 
patno_event_ids = sapply(all_event_ids_p, function(event_id){
  return(paste0(sm$PATNO,'_', event_id ))
})

patno_event_ids=unlist(patno_event_ids)
# select data for the requested patiennts 

combined_bl_log_sel<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids )
combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel$INEXPAGE=='INEXPD',]




correlate_factors_categorical<-function(y){
  ##  obtain correlations of factors with the clinical clusters 
  #' explicitly done for categorical variables
  #' @param name 
  #'
  #'
  #'
  #'
    clust_name=paste0(y, '_clin_clust')
    
    to_clust<-combined_bl_log_sel_pd_to_clust
    to_clust<-combined_bl_log_sel_pd
    
    cl_clusters<-try(get_clinical_clusters_kml(to_clust,y, nbCluster = 5), 
                     silent = TRUE)
    
    if (!class(cl_clusters) == "try-error"){
      cl_clusters_MOFA<-cl_clusters[match( samples_metadata(MOFAobject)$PATNO, names(cl_clusters))]
      names(cl_clusters_MOFA)<-samples_metadata(MOFAobject)$PATNO
      samples_metadata(MOFAobject)[,clust_name ]=cl_clusters_MOFA
      
      
      
      cl_clusters_MOFA
      samples_metadata(MOFAobject)[,clust_name ]
      
      #samples_metadata(MOFAobject)[,clust_name ]=cl_clusters_MOFA
      #correlate_factors_with_covariates 
      sm<-samples_metadata(MOFAobject)
      if (clust_name %in% colnames(sm)){
        covariates=samples_metadata(MOFAobject)[,clust_name ]
        Z <- get_factors(MOFAobject, factors = 1:15, 
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
      }else{
        print('could not calculate clusters')
      }
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



############ SUBSET PATIENTS ONLY ########################################
# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 

#MOFAobjectBL <- subset_groups(MOFAobject, groups = 1)
PD_samples_only<-MOFAobject@samples_metadata$PATNO_EVENT_ID[MOFAobject@samples_metadata$COHORT_DEFINITION=='Parkinson\'s Disease']

HC_samples_only<-MOFAobject@samples_metadata$PATNO_EVENT_ID[MOFAobject@samples_metadata$COHORT_DEFINITION=='Healthy Control']
HC_samples_only


MOFAobjectPD <- subset_samples(MOFAobject, samples=PD_samples_only)

MOFAobjectHC <- subset_samples(MOFAobject, samples=HC_samples_only)

sel_group=4
groups_names(MOFAobject)

## subset group
if (length(groups_names(MOFAobject))>1){
  MOFAobject_sel<-subset_groups(MOFAobject, sel_group)
  MOFAobjectPD_sel<-subset_groups(MOFAobjectPD, sel_group)

  met<-samples_metadata(MOFAobject_V08)
}else{
  MOFAobjectPD_sel = MOFAobjectPD
  MOFAobject_sel = MOFAobject

}





# CORS PATIENTS ONLY #### 
dir.create(paste0(outdir, '/covariates/'))

stats<-apply(MOFAobjectPD_sel@samples_metadata, 2,table )
non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )

covariates_dir<-paste0(outdir, '/covariates/covariate_corelations_')
covars_f_pearson<-paste0(covariates_dir, 'pearson.csv' )
covars_f_pvalue<-paste0(covariates_dir, 'pvalue.csv' )
covars_f_pearson_pd<-paste0(covariates_dir, 'pearson_pd.csv' )
covars_f_pvalue_pd<-paste0(covariates_dir, 'pvalue_pd.csv')

 #names(non_na_vars)
 force_cors=FALSE 
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
covariates=c('SITE', 'Plate', 'NHY', 'NP2PTOT_LOG', 'AGE_SCALED', 'SEX', 
             'Usable.Bases....', 'updrs2_score',  'scopa')
measured_cells<-c('Neutrophils....', 'Lymphocytes....', 'Neutrophil.Score')


#correlate_factors_with_covariates_categ(MOFAobjectPD, covariates=covariates)
sel_factors_pd_np3<-which(cors_all_pd[,c('NP3TOT_LOG' )]>(-log10(0.05)))


### DIFF ####
all_diff<-colnames(sm)[grep('diff', colnames(sm))]
all_diff_variables<-colnames(sm)[grep('diff', colnames(sm))]
all_diff_variables<-colnames(sm)[grep('diff', colnames(sm))]

# variables to cluster 

all_diff_variables=c(all_diff_variables,'NP1RTOT', 'NP2PTOT','NP2PTOT_LOG', 'NP2PTOT_LOG_V10', 'NP3TOT','NP3TOT_LOG', 'updrs3_score', 'updrs2_score',
                     'updrs2_score_LOG', 'updrs3_score_LOG', 
                     'scopa', 'rem', 'upsit', 'moca', 'sft',
                     'abeta', 'sft_V12')

                     # variables to create clusters 
all_diff_variables=c('NP1RTOT', 'NP2PTOT','NP2PTOT_LOG', 'NP2PTOT_LOG_V10', 'NP3TOT','NP3TOT_LOG', 'updrs3_score', 'updrs2_score',
                     'updrs2_score_LOG', 'updrs3_score_LOG', 'updrs3_score_on','updrs3_score_on_LOG',
                     'scopa', 'rem', 'upsit', 'moca', 'sft',
                     'abeta', 'sft_V12')


all_diff_in_cors<-all_diff_variables[all_diff_variables %in% colnames(cors_all_pd)]
all_diff_in_cors<-all_diff_in_cors[!grepl('clust', all_diff_in_cors)]
#round(cors_all_pd[,'moca'], digits=2)
#round(cors_all_pd[,'NP2PTOT_LOG'], digits=2)
#round(cors_pearson_pd[,'NP2PTOT_LOG'], digits=2)

#round(cors_all_pd[,'NP3TOT_LOG'], digits=2)
#round(cors[,'COHORT'], digits=2)
#round(cors_pearson[,'COHORT'], digits=2)

# HERE CHOOSE THE FACTORS THAT ACTUALLY ASSOCIATE with the longterm differences 


























