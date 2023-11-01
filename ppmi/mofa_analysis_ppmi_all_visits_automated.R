
#install.packages('psych')


#install.packages("remotes")
#remotes::install_github("bioFAM/MOFAdata")
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


#as.factor(tolower(samples_metadata(MOFAobject)$SCAU26CT))
#as.factor(samples_metadata(MOFAobject)$SCAU26CT)
MOFAobject@samples_metadata=meta_merged_ord

length(MOFAobject@samples_metadata$PATNO_EVENT_ID)
samples_metadata(MOFAobject)$SCAU26CT<-as.factor(tolower(samples_metadata(MOFAobject)$SCAU26CT))


samples_metadata(MOFAobject)$months<-unlist(EVENT_MAP[samples_metadata(MOFAobject)$EVENT_ID], use.names = FALSE)


### ADD FUTURE CHANGES ####
############################




sm<-samples_metadata(MOFAobject)
res<-df_change[match(sm$PATNO, np3_diff$PATNO),]

samples_metadata(MOFAobject) = cbind(samples_metadata(MOFAobject),res )



##################

library(ggplot2)

#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)

graphics.off()

jpeg(paste0(outdir, 'factor_cor','.jpeg'))
plot_factor_cor(MOFAobject)
dev.off()
group=1

### Change mofa names to gene symbols 
### create a new object with rna names 
MOFAobject_gs<-MOFAobject

ens_ids_full<- features_names(MOFAobject)$RNA
ens_ids<-gsub('\\..*', '', ens_ids_full)

symbols_ids<-get_symbols_vector(ens_ids)
features_names(MOFAobject_gs)$RNA<-symbols_ids
features_names(MOFAobject_gs)$RNA

  
library(ensembldb)
#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)

get_top_cors<-function(MOFAobject, COHORT_NAME='CONCOHORT'){
  cors_both<-get_correlations(MOFAobject, names(non_na_vars))
  cors_top<-cors_both[[1]]
  cors_pearson_top<-cors_both[[2]]
  sel_factors<-abs(cors_pearson_top[,COHORT_NAME])>0.15
  
  round(cors_top[,COHORT_NAME][sel_factors], digits=2)
  round(cors_pearson_top[,COHORT_NAME][sel_factors], digits=2)
}



######### VARIANCE EXPLAINED ###########

group=1


vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]

vars_by_factor_f<-format(vars_by_factor*100, digits=2)
vars_by_factor
write.table(format(vars_by_factor_f,digits = 2)
            ,paste0(outdir,'variance_explained.txt'), quote=FALSE)

p3<-plot_variance_explained(MOFAobject, max_r2=20)+
  theme(axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))
ggsave(paste0(outdir, 'variance_explained','.png'), plot=p3,
       width = 5, height=N_FACTORS/2,
       dpi=100)

MOFAobject@samples_metadata$Outcome

### Data overview
plot_data_overview(MOFAobject)+
  theme(axis.title.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        text  = element_text(size=16))
ggsave(paste0(outdir, 'data_overview.jpeg'), dpi=100, 
       width=4, height=3)




############ SUBSET PATIENTS ONLY ####
# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 

MOFAobjectBL <- subset_groups(MOFAobject, groups = 1)

#
PD_samples_only<-MOFAobject@samples_metadata$PATNO_EVENT_ID[MOFAobject@samples_metadata$COHORT_DEFINITION=='Parkinson\'s Disease']

MOFAobjectPD <- subset_samples(MOFAobject, samples=PD_samples_only)

# CORS PATIENTS ONLY #### 

stats<-apply(MOFAobjectPD@samples_metadata, 2,table )
non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )
cors_both<-get_correlations(MOFAobjectPD, names(non_na_vars))
cors_pearson_pd=cors_both[[2]]; cors_pd=cors_both[[1]]; cors_all_pd=cors_both[[1]]


#### Covariance of factors with metadata 
source('ppmi/mofa_utils.R')


stats<-apply(MOFAobject@samples_metadata, 2,table )
non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )
cors_both<-get_correlations(MOFAobject, names(non_na_vars))
cors_pearson=cors_both[[2]]; cors=cors_both[[1]]; cors_all=cors_both[[1]]


######



if (length(sel_coh)>1){
  sel_factors<-which(cors_all[,'COHORT' ]>-log10(0.05))
  sel_factors_saa<-which(cors_all[,c('CSFSAA' )]>(-log10(0.05)))
  
  
}else{
  sel_factors<-which(cors_all[,c('NP3_TOT' )]>-log10(0.05))
  
}
sel_factors_pd_np3<-which(cors_all_pd[,c('NP3_TOT' )]>(-log10(0.05)))
sel_factors_pd_np3<-which(cors_all_pd[,c('NP3_TOT_LOG' )]>(-log10(0.05)))

sel_factors_pd_np2<-which(cors_all_pd[,c('NP2_TOT_LOG' )]>(-log10(0.05)))






### DIFF ####
sm<-samples_metadata(MOFAobject)
all_diff<-colnames(sm)[grep('diff', colnames(sm))]
all_diff_variables<-colnames(sm)[grep('diff', colnames(sm))]

sel_factors_pd_np3_diff<-which(cors_all_pd[,c('NP3_TOT_diff_V16' )]>(-log10(0.05)))
sel_factors_pd_hi_put_diff<-which(cors_all_pd[,c('hi_putamen_diff_V10' )]>(-log10(0.05)))

# HERE CHOOSE THE FACTORS THAT ACTUALLY ASSOCIATE with the longterm differences 

all_fs_diff<-sapply(all_diff, function(diff_var){
  try(
    
    if (diff_var %in% colnames(cors_all_pd)){
      fs<-print(which(cors_all_pd[,c(diff_var)]>(-log10(0.05))))
        return(fs)
      
    }
  )
})
all_fs_diff
#sel_factors_diff=c('Factor1','Factor4','Factor6', 'Factor8', 'Factor9','Factor11' ,'Factor12'    )

sel_factors_np3<-which(cors_all[,c('NP3_TOT' )]>-log10(0.05))
sel_factors_outcome<-which(cors_all[,'Outcome' ]>-log10(0.05))





### RUN CLUSTERING 
library(DescTools)

#TODO: Adjust mode if cohorts are 1 vs 2





cluster_samples_mofa=TRUE
if (cluster_samples_mofa){
  if (length(sel_coh)>1){
    #for (k_centers_m in c(6)){
    for (k_centers_m in c(3,2)){
      
      clusters <- cluster_samples(MOFAobject, k=k_centers_m, factors=sel_factors)
      clusters_mofa<-clusters
      
      clusters_mofa_34 <- cluster_samples(MOFAobject, k=k_centers_m, factors=c(3,4))
      clusters_mofa_moca <- cluster_samples(MOFAobject, k=k_centers_m, factors=c(14,10,8))
      #clusters_mofa_outcome <- cluster_samples(MOFAobjectBL, k=6, factors=sel_factors_outcome)
      
      clusters_patients <- cluster_samples(MOFAobjectPD, k=k_centers_m, factors=sel_factors)
      
      ### cluster patients only 
      clusters_patients_pd_np3 <- cluster_samples(MOFAobjectPD, k=k_centers_m, factors=sel_factors_pd_np3)
      
    
      
    }
    # clusters <- cluster_samples(MOFAobject, k=3, factors=c(3,4))
    
  }
  MOFAobject@samples_metadata$clusters=clusters$cluster
  
  ss_scores<-c()
  for (k in 3:15){
    clusters_test <- cluster_samples(MOFAobject, k=k, factors=c( 2:14))
    cluster_bt<-clusters_test$betweenss/clusters_test$totss
    print(cluster_bt)
    ss_scores<-append(  ss_scores,cluster_bt)
  }
  plot(ss_scores)
  
  
  
  
  ########### Add some metadata ####
  samples_metadata(MOFAobject)$PATNO_EVENT_ID
  samples_metadata(MOFAobject)$cluster<-factor(clusters_mofa$cluster)
  samples_metadata(MOFAobject)$cluster<-factor(clusters_mofa$cluster)
  
  
}



####
###DO THE CLUSTERS HAVE DIFFERENT NP3? ####
samples_metadata(MOFAobjectPD)$clusters_np3<-factor(clusters_patients_pd_np3$cluster)
met<-samples_metadata(MOFAobjectPD)

y='NP3_TOT'
y='NP3_TOT'
met%>%
  group_by(clusters_np3) %>%
  summarize(across(everything(), mean))
dev.off()
p<-ggplot(met ,aes_string(x='clusters_np3' , fill='clusters_np3'))+
  geom_boxplot(aes_string(y=y, x='clusters_np3'))+
  geom_signif(comparisons=list(c("1", "2")),
                               aes_string(y=y))


p
ggsave(paste0(outdir,'/clustering/', y, '.png'))
## TODO: WILCOX TEST BY GROUP


#### After adding clusters redo calculcations of metadata..? ####




mt_kv<-read.csv(paste0(output_files, 'metadata_key_value.csv'), header = FALSE)
mt_kv$V1<-gsub(' |\'|\"','',mt_kv$V1 )
mt_kv$V2<-gsub(' |\'|\"','',mt_kv$V2 )


ids_to_plot_cor<-colnames(cors_pearson[,colSums(abs(cors_pearson)>0.2)>0L])
ids_to_plot<-which(apply(cors_pearson, 2, sum)>0)

ids_to_plot<-which(apply(cors, 2, sum)>-log10(0.05))

which(cors_pearson>0.1)
ids_to_plot<-which(apply(cors, 2, sum)>0)
ids_to_plot
tot_cor_t=5
ids_to_plot_strict<-which(apply(cors, 2, sum)>-log10(0.0001))
non_na_ids_to_plot<-intersect(names(non_na_vars),names(ids_to_plot) )
non_na_ids_to_plot




format(cors_pearson[, 'CONCOHORT'], digits=1)
format(cbind(cors_pearson[, 'CONCOHORT'],10^-cors_all[,'CONCOHORT']), digits=3) [sel_factors,]

MOFAobject@samples_metadata$DATSCAN_CAUDATE_L

#### All associations that are not NA 


to_remove_covars<-grepl( 'DATE|REC_ID|UPDATE|ORIG_ENTR|INFO', non_na_ids_to_plot)
non_na_ids_to_plot_cleaned<-non_na_ids_to_plot[!to_remove_covars]




#ids_to_plot_strict
#ids_to_plot_strict=c('SEX', 'AGE_AT_VISIT')
graphics.off()
to_remove_regex<-'DATE|REC_ID|UPDATE|ORIG_ENTR|INFO'
to_remove_covars<-grepl( to_remove_regex, names(ids_to_plot_strict))
to_remove_covars
ids_to_plot_strict_cleaned<-ids_to_plot_strict[!to_remove_covars]
to_remove_covars
ids_to_plot_strict_cleaned
ids_to_plot_strict_cleaned
jpeg(paste0(outdir, '/factors_covariates_only_nonzero_strict_pval','.jpeg'),
     width =700+length(ids_to_plot_strict)*20,
     height = 800, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(ids_to_plot_strict_cleaned), 
                                  plot = "log_pval"
                                  
)
dev.off()



names(non_na_vars)[ids_to_plot]
MOFAobject@samples_metadata$SNCA_rs356181

#### correlations between factors 
#cors[,c('stai_state' )]
format(cors_pearson[,c('stai_trait')], digits=2)

#cors[,c('NP1RTOT', 'NP1_TOT','NP2PTOT', 'NP2_TOT','NP3TOT'  ,'NP3_TOT','NP4TOT', 'NP4_TOT' , 'SCAU', 'SCAU_TOT', 'RBD_TOT' )]
#format(cors_pearson[,c('NP1RTOT', 'NP1_TOT','NP2PTOT', 'NP2_TOT','NP3TOT'  ,'NP3_TOT','NP4TOT', 'NP4_TOT' , 'SCAU', 'SCAU_TOT','RBD_TOT' )], digits=2)


selected_covars<-c('COHORT', 'AGE_AT_VISIT', 'SEX', 'NP1TOT', 'NP3TOT', 'NP4TOT', 'SCAU', 'PDSTATE')
samples_metadata(MOFAobject)$Outcome

cors_all[, diff_variables]

sm<-samples_metadata(MOFAobject)
sm$hvlt_immediaterecall
### which other variables relate to np3 tot factors 
## what other variables do np3 factors related with? 
f1_np3<-sel_factors_pd_np3[1]

cors_all[f1_np3,][cors_all[f1_np3,]>2]



#selected_covars

if (length(sel_coh)>1){
  selected_covars2<-c(selected_covars2, 'COHORT')
}

sm$NP3_TOT
selected_covars_img<-c('Disease status','hi_caudate', 'ips_caudate', 'con_putamen', 'con_putamen_V10' )

sm<-samples_metadata(MOFAobject)
sm$asyn
MOFAobjectPD
selected_covars=selected_covars2
MOFAobject_to_plot=MOFAobjectPD
plot_covars_mofa<-function(selected_covars, fname, plot, factors,labels_col=FALSE, height=1000, MOFAobject_to_plot=MOFAobject){
  
  # filter if some do not exist in the colnames of metadata
  #apply(MOFAobjectPD@samples_metadata[,selected_covars3], 2, function(x) {length(which(duplicated(x)))==length(x)-1 })
  # first check if the requested names exist in the metadata 
  selected_covars=selected_covars[selected_covars %in% colnames(MOFAobject_to_plot@samples_metadata) ]
  
  
  sds<-apply(MOFAobject_to_plot@samples_metadata[,selected_covars], 2, sd, na.rm=TRUE)
  sd_na<-c(is.na(sds)|sds==0)
  
  print(selected_covars)
  # then check that the sd is not NA
  selected_covars=selected_covars[ !(sd_na) ]
  
  if (labels_col){
    labels_col<-mt_kv$V2[match(selected_covars,mt_kv$V1)]
    labels_col[is.na(labels_col)]<-selected_covars[is.na(labels_col)]
  }else{
    labels_col=selected_covars
  }
  
  
  
  
  jpeg(paste0(outdir,'/', fname,'.jpeg'), width = 1000+length(selected_covars)*20, height=height, res=
         100)
  P2<-correlate_factors_with_covariates(MOFAobject_to_plot,
                                        covariates =selected_covars , plot = plot,
                                        labels_col=labels_col, 
                                        factors = factors, 
                                        cluster_cols=TRUE)
  
  dev.off()
  
  
}

           #  'DYSKIRAT')
# 'STAID:ANXIETY_TOT'

MOFAobject_gs2<-MOFAobject

graphics.off()

plot="log_pval"
# Plot 1: strict ones we are interested in
# conference poster 


factors=names(sel_factors)
sel_factors
fname<-'factors_covariates_only_nonzero_strict_PD'

plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot,factors=sel_factors,labels_col=FALSE, MOFAobject=MOFAobjectPD )

fname<-'factors_covariates_only_nonzero_strict_PD_np3'
plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot,factors = sel_factors_pd_np3,labels_col=TRUE, MOFAobject=MOFAobjectPD )


fname<-'factors_covariates_only_nonzero_strict_cor_PD_np3'
plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot='r',factors = sel_factors_pd_np3,labels_col=TRUE, MOFAobject=MOFAobjectPD )
graphics.off()











########################## NEW LOGIC GET ONLY FACTORS WITH DIFF IN V16 ##############
############# AFTER WE ADDE THIS TO MOFA
all_diff<-all_diff_variables[all_diff_variables %in% colnames(cors_all)]
all_cors_diff<-cors_all[, all_diff]

sel_factors_diff<-which(rowSums(all_cors_diff)>0)
all_diff_variables_prog<-c(all_diff_variables, 'AGE', 'SEX')

#### Factors related to the longterm change in scale #####
fname<-'factors_covariates_only_nonzero_strict_PD_diff'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot,factors = sel_factors_diff,labels_col=TRUE, MOFAobject=MOFAobjectPD )
fname<-'factors_covariates_only_nonzero_strict_cor_PD_diff'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot='r',factors = sel_factors_diff,labels_col=TRUE, MOFAobject=MOFAobjectPD )




fname<-'factors_covariates_only_nonzero_strict_cor_PD'
plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot='r',factors,labels_col=TRUE, MOFAobject=MOFAobjectPD )


fname<-'factors_covariates_only_nonzero_strict'
plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot,factors,labels_col=TRUE, MOFAobject=MOFAobject )




# Plot 1: some more non motor that we discovered
samples_metadata(MOFAobject)$rigidity
samples_metadata(MOFAobject)$NP2_TOT


fname<-'factors_covariates_only_nonzero_broad_PD'
plot_covars_mofa(selected_covars_broad,fname,plot,c(1:15),labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD  )

fname<-'factors_covariates_only_nonzero_broad_PD'
plot_covars_mofa(selected_covars_broad,fname,plot,c(1:15),labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD  )

fname<-'factors_covariates_only_nonzero_broad_cor_PD'

plot_covars_mofa(selected_covars_broad,fname,plot='r',c(1:15),labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD  )




ind_re<-which(non_na_ids_to_plot %in% c('DYSKIRAT'))
## this is the othet
# = non_na_ids_to_plot[-ind_re]

##jpeg(paste0(outdir, 'factors_covariates_only_nonzero_cor_logpval','.jpeg'), 
#     width = 700+length(selected_covars)*20, height=1100, res=100)

##correlate_factors_with_covariates(MOFAobject,covariates = selected_covars, 
#                                  plot = "log_pval",
                                 # alpha=0.000000001,
#                                  col.lim=c(-0.4, 0.4))
#dev.off()

cors_pearson[,c('DYSKIRAT')]

MOFAobject_nams<-MOFAobject
vars_by_factor/rowSums(vars_by_factor)*100

hist(MOFAobject@samples_metadata[,'DYSKIRAT'])

selected_covars_pearson<-selected_covars[!grepl('con_putamen', selected_covars) ]; length(selected_covars_pearson)

f_to_plot<-names(sel_factors)
ind_to_update<-colnames(MOFAobject_nams@samples_metadata) %in%selected_covars_pearson
colnames(MOFAobject_nams@samples_metadata)[ind_to_update]
MOFAobject_nams@samples_metadata$scopa

### Corelation not log_pval####

fname<-'factors_covariates_img_cor'
plot_covars_mofa(selected_covars=selected_covars_img,fname,plot='r',factors,labels_col=FALSE, MOFAobject=MOFAobject_nams )

fname<-'factors_covariates_img_pval'
plot_covars_mofa(selected_covars=selected_covars_img,fname,plot='log_pval',factors,labels_col=FALSE, MOFAobject=MOFAobject_nams )




### 2. Write the covariates for each factor to files ####
### filter only the ones that are correlated 

#write.csv(covariate_corelations, paste0(outdir, '/covariate_corelations.csv'))
dir.create(paste0(outdir, '/covariates/'))
write.csv(cors_pearson, paste0(outdir, '/covariates/covariate_corelations_pearson.csv'))
for (fx in 1:N_FACTORS){
  sig<-cors[fx,]>1.5
  c1<-cors[fx,][sig]
   c2<-cors_pearson[fx,][sig]
  c3<-format(cbind(c1,c2), digits=2); c3<-c3[order(c3[,1], decreasing = TRUE),]
  write.csv(c3, 
            paste0(outdir, '/covariates/',fx, '.csv'))
}



view='proteomics'; factor=6

vps=length(MOFAobject@dimensions$D)

fps= as.numeric(MOFAobject@dimensions$K)
views<-names(MOFAobject@dimensions$D)



###########################################################
#### Weights ####
##### WRITE ALL weights for each factor in one file 


### Actually get only factors with higher variance in RNA
dir.create(paste0(outdir, 'top_weights/'))

T=0.3

# TODO: save to zip file!
for (i in seq(1,vps)){
  view=views[i]
  
  cluego1<-paste0(outdir, 'top_weights/top_weights_vals_by_view_CLUEGO_', view, '_T_', T, '.txt')
  
  all_weights1<-MOFA2::get_weights(MOFAobject_gs,
                                  views = view, 
                                  as.data.frame =TRUE)  
  # threshold each? 

  all_weights_filt<-all_weights1[abs(all_weights1$value)>T,]
  ens_ids<-gsub('\\..*', '', all_weights_filt$feature)
  write.csv(ens_ids,cluego1,
            row.names = FALSE, quote=FALSE)

  ### write gene symbols here 
  all_weights1<-MOFA2::get_weights(MOFAobject_gs,
                                   views = view, 
                                   as.data.frame =TRUE)  
  
  
  all_weights_filt<-all_weights1[abs(all_weights1$value)>T,]
  write.table(all_weights_filt,paste0(outdir, 'top_weights/top_weights_vals_by_view_', view, '_T_', T, '.txt'), sep = '\t')
  
  # threshold each? 
  
  
  
  }


outdir
high_vars_by_factor<-vars_by_factor>0.1


#### 2. Save highly weighted features #####
for (i in seq(1,vps)){
  for (ii in seq(1,fps)){
    
    
    #### print only the views with high variance in factor  
    view=views[i]
    factor=ii
    print(view, factor)
    all_weights<-MOFA2::get_weights(MOFAobject_gs,views = view, factors=factor, 
                            as.data.frame =TRUE)
    
    ### get the top highly weighted variables - absolute value
    top<-all_weights[order(abs(all_weights$value), decreasing = TRUE),]
    if (high_vars_by_factor[factor, view]){
      write.table(top,paste0(outdir, 'top_weights/top_weights_vals',factor,'_', view,'.txt'), sep = '\t')
      
    }
    
    

  }
  }

  dev.off()
  p1<-plot_variance_explained(MOFAobject, plot_total = T)[[2]]
  p1<-p1+theme(axis.text.x=element_text(size=16), 
               axis.title.y=element_text(size=16), 
               axis.text.y=element_text(size=16))
  plot(p1)
  ggsave(paste0(outdir, 'variance_explained_total','.png'),plot=p1, 
         width = 3, height=3, dpi=300)
#install.packages('psych')

  
  
#### Get all weights and put it in one file
  #1. Collate in a list  - order significant
  #2. Stack the lists - 
  

  
  
  
  
  
### wHICH VARIABLES correlate with which factors 
pos_cors<-cors>0  # which are sig. corelated with any factors 
n_factors_pos=1
positive_cors<-cors[,colSums(pos_cors)>n_factors_pos]


### THRESHOLD: IMPORTANT
### significance-which variables to print--> log10(pval)=4 - sig is anything above 1.3
x_cor_t=2

i=111
#### Factor plots ####

positive_cors_to_plot<-positive_cors[,!grepl(tolower(to_remove_regex),tolower(colnames(positive_cors))) ]
positive_cors_to_plot
grepl(positive_cors,names(positive_cors_to_plot))
#install.packages('forcats')
select_factors_manually=TRUE
for (i in 1:dim(positive_cors_to_plot)[2]){
      #' fix find 2d plots 
      #' i= index of the clinical variable 
    
      names<-colnames(positive_cors_to_plot)
      x_cors<-positive_cors_to_plot[,i]
      print(names[i])
      ### does the variable relate to two factors? 
      pos_factors<-names(which(x_cors > 0))
      pos_factors<-names(which(x_cors > x_cor_t))
      pos_factors
      
      
      # Order by 
      pos_factors<-pos_factors[order(x_cors[pos_factors], decreasing = TRUE)]
      print(paste(i, pos_factors))
      
      if (select_factors_manually){
        pos_factors=c(1,4)
      }
      if (length(pos_factors)){
      
              #TODO: print also other factors combinations
              #combn(pos_factors,2)
              
              if (length(pos_factors)>1){fs=c(pos_factors[1],pos_factors[2])}else{fs=pos_factors[1]}
              if  (length(pos_factors)>1){
              
                
                      color_by<-names[i] ## OR change: COLOR BY GROUP? color_by='group'
                      print(color_by)
                      factor_cors<-paste0(format(x_cors[pos_factors], digits=2), collapse=',')
            
                      
                      pf<-plot_factors(MOFAobject, 
                                   factors = fs, 
                                   #shape_by=color_by,
                                   color_by = color_by,
                                     show_missing = FALSE 
                      )
                      
                      pf=pf+labs(caption=paste0('log10pval = ',factor_cors))
                      pf
                      fss<-paste(fs,sep='_',collapse='-')
                      dir.create(file.path(paste0(outdir,'/factor_plots/2D/')), showWarnings = FALSE)
                      
                      FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by, x_cor_t,'.png')
                      
                      
                      ggsave(FNAME,plot=pf, width = 4, height=4, dpi=100)
              
              }
              ### you can also add a variable if it is also related to these two factors? 
              #shape_by='NHY'
              #shape_by='AGE_AT_VISIT'
              #fs
              #pf=plot_factors(MOFAobject, 
             #              factors = fs, 
            #               color_by=color_by,
            #                shape_by= shape_by,
            #               show_missing = FALSE
            #  )
            #  pf=pf+labs(caption=paste0('log10pval = ',factor_cors))
            #  
            #  FNAME<-paste0(outdir,'/factor_plots/group/2D/', 'plot_factors_variate_2D',fss,'_',color_by,'_',shape_by, x_cor_t,'.png')
            #    ggsave(FNAME,plot=pf, width = 4, height=4, dpi=100)
      }
}
dev.off()


######## Specific correlations #########
library(tidyverse)
factor_cors<-paste0(format(x_cors[factors_to_plot], digits=2), collapse=',')
factor_cors
  
MOFAobject@samples_metadata$td_pigd
factors_to_plot<-c(3,4)
color_by='td_pigd'
color_by='NP3GAIT'

pf<-plot_factor(MOFAobject_gs, 
                 factors = c(3,4), 
                # shape_by=color_by,
                 color_by = color_by,
                 show_missing = TRUE 
)
pf
pf=pf+labs(caption=paste0('log10pval = ',factor_cors))
fss<-paste(fs,sep='_',collapse='-')
dir.create(file.path(paste0(outdir,'/factor_plots/2D/')), showWarnings = FALSE)

FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by, x_cor_t,'.png')



MOFAobject@samples_metadata$NP3_TOT







########## Scatter plots - FACTORS to variables ################
# here find the view for which the variability of the factor maximum
plot_data_scatter_by_factor<-function(factor, color_by,MOFAobject_gs=MOFAobject){
  #'
  #' @param 
  #'
  
  top_view<-which.max(vars_by_factor[factor,])
  plot_data_scatter(MOFAobject_gs, view = top_view,factor = factor,  features = 15,sign = "negative",color_by = color_by) + 
    labs(y=color_by)
  ggsave(paste0(outdir, '/scatter_plots/',  views[top_view], '_', factor, '_',  color_by, '.png' ), height=8, width=12, dpi=100)
  
  
}

color_by='NP3_TOT';
dir.create(paste0(outdir, '/scatter_plots/'))
sel_factors_np3
sapply(sel_factors, plot_data_scatter_by_factor, color_by=color_by, MOFAobject_gs=MOFAobjectPD)




### Scatter plot top features NP3
# create function to give top n
top_mirs<-get_weights(MOFAobject_gs, view='miRNA', as.data.frame = TRUE, factor=3)

view='miRNA'
top_mirs<-get_weights(MOFAobject_gs, view=view, as.data.frame = TRUE, factor=3)

top_mirs_10<-top_mirs[order(abs(top_mirs$value), decreasing = TRUE), ][1:20,]
top_mirs_10$feature
mirdata<-MOFAobject_gs@data[view]$miRNA[[1]]


mirdata_sel<-mirdata[rownames(mirdata)%in%top_mirs_10$feature,]
rownames(mirdata_sel)
mirdata_sel_t<-t(mirdata_sel)

MOFAobject@samples_metadata$PDSTATE

## Loop cl 
CL='NP3_TOT'
mirdata_sel_t<-data.frame(mirdata_sel_t, 
                          COHORT=MOFAobject@samples_metadata$COHORT, 
                          CL=MOFAobject@samples_metadata[, CL], 
                          PDSTATE=MOFAobject@samples_metadata[, 'PDSTATE'])


MOFAobject@samples_metadata$COHORT_DEFINITION
colnames(mirdata_sel_t$hsa.miR.193b.3p)
mirdata_sel_t$hsa.miR.193b.3p
mirdata_sel_t$PATNO<-rownames(mirdata_sel_t)
mirdata_sel_t<-mirdata_sel_t[mirdata_sel_t$COHORT==1,]
mirdata_sel_t$PDSTATE
mirdata_sel_t<-mirdata_sel_t[!mirdata_sel_t$PDSTATE %in% c('ON')  ,]



library('ggpmisc')

mirdata_melt<-melt(mirdata_sel_t, id=c('PATNO', 'COHORT', 'CL', 'PDSTATE'))
colnames(mirdata_melt)

ggplot(mirdata_melt, aes_string(x='value', y='CL'))+
  geom_point()+
  geom_smooth(method=lm)+
  
  stat_fit_glance(method = 'lm',
                  method.args = list(),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.35, size = 3)

  facet_wrap(~variable, scales='free_x')


ggsave(paste0(outdir, '/trajectories/clinical/',view,CL, '.jpeg'  ), width=10,height=10, dpi=300)


type(MOFAobject@samples_metadata$STAIAD3)

## todo extract row and column for which positive cors are TRUE 
#positive_cors_1<-pos_cors>0
#which(positive_cors_1)

#for (i in 1:length(positive_cors_1)){
#  
#}

color_by='COHORT'
MOFAobject@dimensions$K

plot_factors(MOFAobject, 
             factors =  c(1,4), 
             color_by = color_by,
             
             show_missing = FALSE
)

color_by='COHORT'

cors_both




### Age, gender , stage does not discriminate factors
#Conclusion here:# factor 1 correlates with grade
for (ii in seq(1,fps)){
  ### Plot factors against a clinical variable 
  x_cor_t=4
  cors_sig<-names(which(positive_cors[ii,]>x_cor_t))
  ln_cs<-length(cors_sig)
  if (ln_cs>0){
    for (iii in seq(1:ln_cs)){
         color_by=cors_sig[iii]
          p<-plot_factor(MOFAobject, 
                         factors = ii, 
                         color_by =color_by,
                         add_violin = TRUE,
                         dodge = TRUE,
                         show_missing = FALSE
                        
          )
        
        
        FNAME<-paste0(outdir,'/factor_plots/', 'plot_factors_variate_1D_',ii, '_',color_by,'_cor_', x_cor_t, '.png')
        
        
        ggsave(FNAME, width = 4, height=4, dpi=100)
        
        
    }
  }
  
  
}
cors_pearson[,'COHORT']
f1<-get_factors(MOFAobject,factors=1)
f1<-as.data.frame(get_factors(MOFAobject,factors=1)['group1'])

var_name='STAIAD22'
yvar<-MOFAobject@samples_metadata[var_name]
yvar_name<-'yvar'
length(yvar)
f1[yvar_name]<-yvar
f1[,yvar_name]=as.numeric(f1[,yvar_name])
ggplot(f1, aes_string(x='Factor1', y=yvar_name) )+ geom_point()


# Factor 2 associates with proteomic Subtype 

color_by<-'NHY';fs<-c(4,8)
color_by<-'NP3TOT';fs<-c(1,4)


plot_factor(MOFAobject, 
            factors = fs, 
            color_by = color_by,
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)
fss<-paste(fs,sep='_',collapse='-')
ggsave(paste0(outdir, 'plot_factor_variate_violin',fss,color_by,'.png'), width = 4, height=4, dpi=100)

#### plot 2 factors 
##### TODO: Plot only significant covariates  here

color_by<-'NHY'
fs<-c(4,8)

color_by<-'NHY';fs<-c(1,8)
color_by<-'NP3TOT';fs<-c(1,8)

plot_factors(MOFAobject, 
            factors = fs, 
            color_by = color_by,
            show_missing = FALSE
)
fss<-paste(fs,sep='_',collapse='-')
FNAME<-paste0(outdir, 'plot_factors_variate_2D',fss,color_by,'.png')
FNAME
color_by
ggsave(FNAME, width = 4, height=4, dpi=100)



##### plot weights ####
library(grid)
library(gridExtra)
v_set=c()

view='miRNA'

fps=MOFAobject@dimensions$K
vps
fps
seq(1,fps)
seq(1,vps)
factor=3; view='RNA'
plot_top_weights(MOFAobject,
                 view = view,
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave(paste0(outdir, 'top_weights_',factor, view,'_','.png'), width =3 , height=4, dpi=300)
dir.create(paste0(outdir, 'top_weights/'))


graphics.off()

#####################
#### 3. Save heatmaps and top weighted feaqtures ####
dir.create(paste0(outdir, '/heatmap/'))
views[i]
i
ii=4
i=2
views[i]


data.frame()


#### Load features 
fs_metadata<-read.csv(paste0(data_dir,'/ppmi/ppmi_data/features_metadata_genes.csv'))
colnames(fs_metadata)<-c('f', 'feature_id', 'known')
colnames(fs_metadata)
fs_metadata$view='RNA'
fs_met_to_merge<-fs_metadata[,c('feature_id', 'view', 'known'  )]

source(paste0(script_dir, 'ppmi/mofa_my_utils.R'))
p_ws<-plot_top_weights2(MOFAobject_gs,
                       view = 'RNA',
                       factor = 14,
                       nfeatures = 15,     # Top number of features to highlight
                       scale = F
)
p_ws
# known genes 
for (i in seq(1,vps)){
  for (ii in seq(1,fps)){
   
    nFeatures=12
    
    ### oNLY SAVE THE ones with high variance
    if (high_vars_by_factor[ii, i]){
      print(c(i,ii))
      cols <- c( "red", 'red')
          p_ws<-plot_top_weights(MOFAobject_gs,
                           view = i,
                           factor = c(ii),
                           nfeatures = nFeatures,     # Top number of features to highlight
                           scale = F
          )
          graphics.off()
          if (views[i]=='RNA'){
            p_ws<-plot_top_weights2(MOFAobject_gs,
                                    view = 'RNA',
                                    factor = ii,
                                    nfeatures = nFeatures,     # Top number of features to highlight
                                    scale = F
            )

          
            print('rna')
            p_ws<-p_ws+ theme(axis.text.y = element_text(face = "italic"))
          }
          
          if (views[i]=='metabolites'){
           p_ws<-p_ws+ theme(axis.text.y = element_text(size=7))
          }
          
          #scale_colour_(values=cols) 
       
            
          ggsave(paste0(outdir, 'top_weights/top_weights_','f_', ii,'_',views[i],'.png'), 
                 plot=p_ws, 
                 width = 3, height=nFeatures/4, dpi=300)
          
          plot_weights(MOFAobject_gs, 
                       view = views[i], 
                       factor = ii, 
                       nfeatures = 30
          )
         
            
          ggsave(paste0(outdir, 'top_weights/all_weights_','f_', ii,'_',views[i],'.png'),
                 width = 3, height=nFeatures/4, dpi=300)
      
          

         
          }
  }
}




    # top weights
    # concat all 
    
    
 


#### 3. Save heatmaps and top weighted feaqtures ####


vps
fps
# rename because value is too long in the legend
MOFAobject@samples_metadata$CONCOHORT_DEFINITION[MOFAobject@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'
MOFAobject_gs@samples_metadata$CONCOHORT_DEFINITION[MOFAobject_gs@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'

graphics.off()
for (i in seq(1,vps)){
  for (ii in seq(1,fps)){
    print(paste('Modality', i, 'factor', ii))
    cluster_rows=TRUE;cluster_cols=TRUE
    
    
    
    ###### Heatmaps 
    nfs=40
    #jpeg(paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs, '_cr_',cluster_rows, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
    # Plot heatmaps for each factor only for miRNA 
    
    var_captured<-round(vars_by_factor[ii,i], digits=2)
    main_t<-paste0('Factor ', ii, ', Variance = ',var_captured, '%')
    #log10(0.005) = 2.3
    #cor_T<- -log10(0.005); cor_p_T<-0.15
    modality=names(MOFAobject@dimensions$D)[i]
    if (names(MOFAobject@dimensions$D)[i]=='proteomics'){modality=paste(TISSUE, modality )  }
    main_t<-paste0('Factor ', ii, ', Mod ',modality, ', Variance = ',var_captured, '%')
    
    
    
    ns<-dim(MOFAobject@samples_metadata)[1]
    if (run_mofa_complete){
      cor_T<-1.5; cor_p_T<-0.1
      
    }else{
      cor_T<-2; cor_p_T<-0.1
      
    }
    
    abs(cors_pearson)>0.15
    
    rel_cors<-cors[ii,][cors[ii,]>cor_T &  abs(cors_pearson[ii,])>cor_p_T ]

    # sig holds the names only 
    cors_sig=names(rel_cors); cors_sig
    FT=0
    if (length(cors_sig)==0){
      cors_sig=c()
      
    } else if (length(cors_sig)>10){
      FT=10
      # rel_cors_ordered<-rel_cors[order(-rel_cors)][1:7]
      rel_cors_ordered<-rel_cors[order(-rel_cors)]
       rel_cors_ordered<-rel_cors[order(-rel_cors)][1:FT]
      #rel_cors_ordered<-rel_cors[order(-rel_cors)]

      cors_sig<-names(rel_cors_ordered)
    }
    cors_sig
    exclude_vars= c('LAST_UPDATE_M4', 'INFODT_M4', 'NTEXAMTM', 'REC_ID_moca', 'REC_ID_st')
    #'OFFEXAMTM', 
     #               'OFFEXAMDT', 'OFFPDMEDT', 'INFO')
    cors_sig<-cors_sig[!(cors_sig %in% exclude_vars)]; cors_sig
    cors_sig<-cors_sig[!grepl( 'LAST_UPDATE|INFO_DT|TM|DT|ORIG_ENTRY|DATE|PAG_', cors_sig)]
    
    plot_heatmap_flag=TRUE
    #MOFAobject_gs@samples_metadata[cors_sig][is.na(MOFAobject_gs@samples_metadata[cors_sig])]<-10^-6v
    MOFAobject_gs@samples_metadata[,cors_sig]
    MOFAobject_gs@samples_metadata[,cors_sig]
    
    #is.na(MOFAobject_gs@samples_metadata[,cors_sig])

    ### if the col contains only NA
    
    #which(cors_sig_non_na=='PDSTATE')
    #cors_sig_non_na=cors_sig_non_na[-3]
    if (length(cors_sig)>1){
      cors_sig_non_na<-names(which( !apply(is.na(MOFAobject_gs@samples_metadata[,cors_sig]),2,any )))
      
    }else{
      cors_sig_non_na=cors_sig 
    }
    cors_sig_non_na
    if( length(cors_sig_non_na)==0){
      cors_sig_non_na=c()
    }
    denoise=FALSE
    
    #cors_sig_non_na=cors_sig
    groups='all';groups=2; 
    groups=1
    #hname<-paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, res, '_cor_', cor_T, 'FT_', FT, '.jpeg')
    hname<-paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, '_cor_', cor_T, 'FT_', 
                  FT, 'den_', denoise, groups, '.jpeg')

    #View(MOFAobject_gs@samples_metadata[cors_sig_non_na])
    p<-plot_data_heatmap(MOFAobject_gs, 
                         view = views[i], 
                         factor =  ii,  
                         features = nfs,
                         groups = groups, 
                         denoise = denoise,
                         cluster_rows = cluster_rows, 
                         cluster_cols = cluster_cols,
                         show_rownames = TRUE, show_colnames = TRUE,
                         scale = "row",
                         annotation_samples = cors_sig_non_na,
                         main=main_t
                         
                         
    )
    #ggsave(hname, plot=p,height=nfs/2, width=(ns+as.numeric(length(cors_sig_non_na) )) )
    if (run_mofa_complete){
      width=ifelse( length(cors_sig_non_na)> 0,ns/10+6,ns/10+4)
      
    }else{
      width=ifelse( length(cors_sig_non_na)> 0,ns/80+6,ns/80+4)
      
    }
    
    ggsave(hname, plot=p,height=nfs/5+2, width=width, dpi=250) 
    
    
  }
  # top weights
  # concat all 
  
  
  
}

views[3]

p<-plot_data_heatmap(MOFAobject_gs, 
                     view = views[3], 
                     factor =  ii,  
                     features = nfs,
                     denoise = FALSE,
                     cluster_rows = cluster_rows, 
                     cluster_cols = cluster_cols,
                     show_rownames = TRUE, show_colnames = TRUE,
                     scale = "row",
                     annotation_samples = cors_sig_non_na,
                     main=main_t
                     
                     
)



# plot heatmaps by view and 
n_groups=length(MOFAobject@data_options$groups)
if (length(MOFAobject@data_options$groups)>1){
  
        
        library(cowplot)
        library(ComplexHeatmap)
        nfs=40
        
        breaksList = seq(-3, 3, by = 0.01)
        groups='all'
        
        groups=c(1,3)
        view=c(3)
        factor=8
        image_path<-paste0(outdir, 'heatmap_',ii,'_',view, 'nfs_', factor, 'gr_', groups, '.jpeg')
        
        # set width of plot based on number of samples retrieved
        n_samples<-get_factors(MOFAobject, factors = 1, groups=groups)
        ns<-length(unlist(n_samples))
        
        jpeg(image_path, height=30*nfs, width=20*ns)
        
        
        p1<-plot_data_heatmap(MOFAobject, 
                          view = view,
                          factor = factor,  
                          features = nfs,
                          groups=groups,
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          show_rownames = TRUE, show_colnames = TRUE,
                          scale = "row"
                          
        )
        p1
        dev.off()

}




samples_order<-MOFAobject@samples_metadata$PATNO_EVENT_ID.x
# ORDER patients by their groups
samples_order<-samples_order[ with(MOFAobject@samples_metadata, order(EVENT_ID, PATNO))]

if (n_groups>1){
  p2<-plot_data_heatmap(MOFAobject, 
                        view = "proteomics",
                        factor = 1,  
                        features = nfs,
                        groups='all',
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        show_rownames = TRUE, show_colnames = TRUE,
                        scale = "row", 
                        max.value = 3, 
                        width=10
  )
  p2
}




library(gridExtra)
#grid.arrange(arrangeGrob(grobs=list(p1, p2), nrow = 1, top="Main Title"))
#do.call('grid.arrange', c(list(p1,p2)) )

#dev.off()


plot_data_heatmap(MOFAobject, 
                  view = "miRNA",
                  factor = 2,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
#
#plot_data_heatmap(MOFAobject, 
#                  view = "proteomics",
#                  factor = 3,  
#                  features = 100,
#                  cluster_rows = FALSE, cluster_cols = TRUE,
#                  show_rownames = TRUE, show_colnames = FALSE,
#                  scale = "row"
#)


### So multi omics factors are more related to Stage than to subtype!!
##How?? plot on the factors and color by stage

p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "NP3BRADY",
                  shape_by = "NP3BRADY",
                  dot_size = 2.5,
                  show_missing = T
)

#p <- p + 
#  geom_hline(yintercept=-1, linetype="dashed") +
#  geom_vline(xintercept=(-0.5), linetype="dashed")
print(p)
ggsave(paste0(outdir,'factor_plot','.png'), width = 4, height=4, dpi=120)





#### Now make predictions
### 



####
##### 5. GENE SET ENRICHMENT! #####
## AND reactome gsa enrichment!!


###### GSEA 

#BiocManager::install('AnnotationHub')

#source('enrichment.R')
  
#library(AnnotationHub)







##### Make predictions 
## Here we select the top factors associated with the clinical variable of interest. 
## then we can choose the top variables of those factors 

## Prediction of clinical subgroups 
## Predict the EORTC.risk

#install.packages('caret')
#install.packages('tibble')
library('tibble')
library('caret')
packages <- c('dplyr', 'ggplot2', 'caret', 'party')
invisible(lapply(packages, library, character.only = TRUE))

suppressPackageStartupMessages(library(randomForest))

# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=1:15)[[1]])
df <- as.data.frame(get_factors(MOFAobject, factors=sel_factors)[[1]])





# Do predictions with the factors 
# Train the model for eortc.risk

# Use the pricnipal components 
# Train-validation split #
# Cross-Validation ####


df$y <- as.factor(MOFAobject@samples_metadata$NHY)

df$y<- as.factor(MOFAobject@samples_metadata$COHORT)
df_age <- cbind(df,MOFAobject@samples_metadata[, c('AGE_SCALED', 'SEX')])



## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 5 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=5, repeats=3)
train_index <- createDataPartition(y=df$y, p=0.7, list=FALSE)

val_folds<-createFolds(df$y, k = 10, list = TRUE, returnTrain = TRUE)
val_folds#


## TODO: issues : there is class imbalance 
# TODO: issues 
res<-sapply(val_folds, cross_val_score, df=df)
res_age<-sapply(val_folds, cross_val_score, df=df_age)











#########





model.COHORT <- randomForest(COHORT ~ ., data=training_set, ntree=10)
MOFAobject@samples_metadata$COHORT.pred <- stats::predict(model.COHORT, df)

# Assess performance 
predicted<-as.factor(MOFAobject@samples_metadata$COHORT.pred)
actual = MOFAobject@samples_metadata$COHORT
confusion_mat = as.matrix(table(actual, predicted )) 




print(confusion_mat)


N_FACTORS
## Show "importance" of variables: higher value mean more important:
round(importance(model.COHORT), 2)
#install.packages('randomForest')

# Prepare data
# Predict EORTC.risk with factor 1,2 only!
high_weights=get_weights()
df <- as.data.frame(get_factors(MOFAobject, factors=c(2,3,4,6))[[1]])
df_genes<-as.data.frame(get_weights(MOFAobject, factors=c(2,3,4,6)) )

# Train the model for IGHV
y_predict='CONCOHORT_DEFINITION'
df$y <- as.factor(MOFAobject@samples_metadata[,y_predict])
model.y <- randomForest(y ~ .,data= df, ntree=10)
df$y <- NULL # important 


# Do predictions
MOFAobject@samples_metadata$y.pred <- stats::predict(model.y, df)
MOFAobject@samples_metadata$y.pred

# Assess performance 
predicted<-MOFAobject@samples_metadata$y.pred
actual <-as.factor(MOFAobject@samples_metadata[,y_predict])
confusion_mat = as.matrix(table(actual, predicted )) 

print(confusion_mat)
View(confusion_mat)
round(importance(model.y), 2)


#install.packages('GGally')
library('GGally')
### Plot predictions
p <- plot_factors(MOFAobject, 
                  factors = c(2,3,4,6), 
                  color_by = "CONCOHORT_DEFINITION",
                  shape_by = "CONCOHORT_DEFINITION",
                  dot_size = 2.5,
                  show_missing = T
)


show(p)
dev.off()
cbind(MOFAobject@samples_metadata$sample,MOFAobject@samples_metadata$COHORT_DEFINITION)
MOFAobject@samples_metadata[MOFAobject@samples_metadata$PATNO=='3156',]




