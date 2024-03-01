

#install.packages('DescTools')
library(DescTools)
library('factoextra')
# Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model <- load_model(file)

# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 
# cluster samples here 


PD_samples_only<-MOFAobject@samples_metadata$PATNO_EVENT_ID[MOFAobject@samples_metadata$COHORT_DEFINITION=='Parkinson\'s Disease']

MOFAobjectPD <- MOFA2::subset_samples(MOFAobject, samples=PD_samples_only)


get_factors_for_metric<-function(diff_var){

  # get associated factors, remove the ones related to confounding
  fact <- which(all_fs_diff[,diff_var])

        if (remove_cell_factors){
          fact<-fact[!(fact %in% fact_neutro_pd)]

        }
       
        return(fact)
        }


#### Create the MOFA clusters with the same K ####
k_centers_m=3
remove_cell_factors = FALSE
#diff_var=y;  # diff_var='NP2PTOT_diff_V16'
rescale_option=TRUE
#diff_var='NP2PTOT_LOG'

DIFF_VAR='moca'
DIFF_VAR='NP2PTOT_LOG'


sel_group_cors = FALSE # Where to get corelations from
if (sel_group_cors){
  MOFAobjectPD_cors<-subset_groups(MOFAobjectPD, sel_group_cors)
}else{
  MOFAobjectPD_cors<-MOFAobjectPD
}


cors_both_clinical<-get_correlations(MOFAobjectPD_cors, all_diff_in_cors)
cors_pearson_pd_clinical = as.data.frame(cors_both_clinical[[2]]);  cors_all_pd_clinical = as.data.frame(cors_both_clinical[[1]])
## choose if the clu
all_fs_diff_all_time<-as.data.frame(cors_all_pd_clinical[,all_diff_in_cors]>(-log10(0.05)))
all_fs_diff = all_fs_diff_all_time
all_fs_diff[, c('NP2PTOT_LOG', 'updrs2_score_LOG_diff_V12', 'moca')]

### Select group for plotting
sel_group=4
if (length(groups_names(MOFAobject))>1){
  MOFAobject_sel<-subset_groups(MOFAobject, sel_group)
  MOFAobjectPD_sel<-subset_groups(MOFAobjectPD, sel_group)

  met<-samples_metadata(MOFAobject_V08)
}else{
  MOFAobjectPD_sel = MOFAobjectPD
  MOFAobject_sel = MOFAobject
  met<-samples_metadata(MOFAobject)

}
cors<-as.data.frame(cors)
# Remove factors that have to do with site, plate, 




fact_neutro_pd<-c(  (which(cors$`Lymphocytes....`>-log10(0.05))), which(cors$`Neutrophils....`>-log10(0.05)),  which(cors$Plate>-log10(0.05)), 
                      which(cors$RIN>-log10(0.05)), which(cors$Usable_bases_SCALE>-log10(0.05)))
                      #which(cors$`Uniquely.mapped....`>-log10(0.05))
                        

    


all_clusts_mofa <- sapply(colnames(all_fs_diff),function(diff_var){
  #'
  #' @param
  #' 
  #' 

        fact <- get_factors_for_metric(diff_var)
    


        #print(paste(diff_var,paste(fact, collapse=', ')))
        xname = paste0(diff_var, '_clust')
      if (length(fact) > 0){
        set.seed(42)
        set.seed(60)
        clusters_x=cluster_by_mofa_factors(MOFAobject=MOFAobjectPD_sel, centers=k_centers_m,
         factors=fact, rescale=rescale_option) # Cluster using kmeans
        clusters_x
       # if (MOFAobject==1){
        #  names(clusters_x)<-gsub('\\_.*', '', names(clusters_x))

       # }

        
        clust_ps<-clusters_x;clust_ps
        MOFAobjectPD_sel@samples_metadata[, xname] = as.factor(clust_ps)
        # TODO: save labels

        return(clust_ps)

}}
)



all_clusts_mofa$NP2PTOT_LOG
all_clusts_mofa_true<-all_clusts_mofa[!as.logical(lapply(all_clusts_mofa, is.null))]

all_clusts_mofa_true_t<-do.call(cbind,all_clusts_mofa_true)

dir.create(paste0(outdir,'/clustering/'))
all_clusts_file<-paste0(outdir,'/clustering/all_clusts_mofa.csv')

write.csv(all_clusts_mofa_true_t,all_clusts_file, row.names=TRUE,sep=','  )
  clusters_ids<-all_clusts_mofa[['NP2PTOT_LOG']]





for (diff_var in names(all_clusts_mofa)){
    MOFAobjectPD_sel@samples_metadata[,paste0(diff_var, '_clust')]<-all_clusts_mofa[[diff_var]];#all_clusts_mofa[[diff_var]]
    sm<-MOFAobject_sel@samples_metadata
    clusters_ids<-all_clusts_mofa[[diff_var]]
    MOFAobject_sel@samples_metadata[,paste0(diff_var, '_clust')]<-clusters_ids[match(sm$PATNO_EVENT_ID,rownames(clusters_ids ) )]
    MOFAobject_sel@samples_metadata[(sm$INEXPAGE %in% c('INEXHC')),paste0(diff_var, '_clust')]<-'HC'
}





### Boxplots by cluster 
## Can produce for multiple metrics 

diff_variables_to_p=c('NP2PTOT_LOG', 'NP2PTOT','scopa','updrs3_score', 
                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
                      'upsit', 'VLTFRUIT', 'sft', 
                      'stai_state', 'stai_trait', 
                      'AGE' )
                      


#other_metrics<-t(cors_all_pd[all_fs_diff[,y],])

diff_variables_to_p=c('NP2PTOT', 'Neutrophil.Lymphocyte', 'AGE_SCALED', 'scopa', 
         'Neutrophils....', 'Lymphocytes....' ,   'sft', 'hvlt_immediaterecall',  # current 
         'updrs3_score_LOG_V12',
          'hvlt_immediaterecall_V12'

         )

