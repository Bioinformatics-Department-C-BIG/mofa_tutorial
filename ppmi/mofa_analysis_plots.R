


# get correlations for selected vars 

correlate_factors_with_covariates
cors_pd<-correlate_factors_with_covariates(MOFAobjectPD_sel, c('NP2PTOT_LOG'),return_data=TRUE , alpha=1  )



cors_both_clinical<-get_correlations(MOFAobjectPD, all_diff_in_cors)
cors_pearson_pd_clinical = as.data.frame(cors_both_clinical[[2]]);  cors_all_pd_clinical = as.data.frame(cors_both_clinical[[1]])
## choose if the clu


all_fs_diff_all_time<-as.data.frame(cors_all_pd_clinical[,all_diff_in_cors]>(-log10(0.05)))
all_fs_diff = all_fs_diff_all_time
all_fs_diff[, 'NP2PTOT_LOG']

#### RUN CLUSTERING ####
library(DescTools)

#TODO: Adjust mode if cohorts are 1 vs 2

if (add_clinical_clusters){
  cors_kruskal_log<-(-log10(cors_kruskal))
  th<-(-log10(0.05))
  cors_kruskal_log[cors_kruskal_log < th ] <-0
  
  pheatmap(cors_kruskal)
}



######### Cluster patients with clinical trajectory factors ####
### get metrics to test how many clusters to use ####

#BiocManager::install('M3C')
#library(M3C)
library('cluster')

y='abeta'
y='NP2PTOT'
y='updrs2_score_LOG'
y='NP2PTOT_LOG'

factors_to_clust<-which(all_fs_diff[ ,y]);print(factors_to_clust)
sm<-samples_metadata(MOFAobjectPD)
all_clusts<-colnames(sm)[grep('clust',colnames(sm))]
cluster_samples_mofa=TRUE
set.seed(123)
set.seed(1239)

Z <- get_factors(MOFAobjectPD, factors = c(factors_to_clust))[[1]];
as.data.frame(Z)
#Z_scaled<-apply(as.data.frame(Z), 2, scale);

Z_scaled<-Z
graphics.off()
fviz_nbclust(Z_scaled, kmeans, method = "wss",  k.max = 15 )#
  #geom_vline(xintercept = 3, linetype = 2)
fviz_nbclust(Z_scaled, kmeans, method = "silhouette", k.max = 15 )

gap_stat <- clusGap(Z_scaled , FUN = kmeans, nstart = 25,
                    K.max =10 , B = 10)
print(gap_stat, method = "globalmax")
fviz_gap_stat(gap_stat)



#### Color the mofa clusters on the factors plot

########### 
all_fs_diff<-as.data.frame(all_fs_diff)
all_fs_diff[,y]



### Gap statistic
library(cluster)
#set.seed(123)
# Compute gap statistic for kmeans
# we used B = 10 for demo. Recommended value is ~500
MOFAobject

if (cluster_samples_mofa){
  if (length(sel_coh)>1){
    #for (k_centers_m in c(6)){
    for (k_centers_m_try in c(3,2)){
      clusters <- cluster_samples(MOFAobject, k=k_centers_m_try, factors=c(1:N_FACTORS))
    }
    # clusters <- cluster_samples(MOFAobject, k=3, factors=c(3,4))
  }
  MOFAobject@samples_metadata$clusters=clusters$cluster
  
  ss_scores<-c()
  
  Z <- get_factors(MOFAobjectPD, factors = c(1,2))
  
  
  
  
  for (k in 3:15){
    clusters_test <- cluster_samples(MOFAobjectPD, k=k, factors=c( 1,12))
    cluster_bt<-clusters_test$betweenss/clusters_test$totss
    print(cluster_bt)
    ss_scores<-append(  ss_scores,cluster_bt)
  }
  plot(ss_scores)
  
  
  
  
  ########### Add some metadata ####
  samples_metadata(MOFAobject)$PATNO_EVENT_ID
  
  
}



####
###DO THE CLUSTERS HAVE DIFFERENT NP3? ####


met<-samples_metadata(MOFAobjectPD_sel)


#boxplot_by_cluster(met, y='PD_MED_USE', clust_name=clust_name )




samples_metadata(MOFAobject_sel)$B.Cells
# Object to hold renamed features - add to metadata..? 
MOFAobject_gs<-MOFAobject_sel
ens_ids_full<- features_names(MOFAobject)$RNA
ens_ids<-gsub('\\..*', '', ens_ids_full)
features_names(MOFAobject_gs)$RNA<-get_symbols_vector(ens_ids)

features_names(MOFAobject_gs)$RNA

MOFAobject_gs2<-MOFAobject

graphics.off()

plot="log_pval"
# Plot 1: strict ones we are interested in
# conference poster 



MOFAobject@samples_metadata$tau_asyn



colnames(cors_all_pd)[grep('tau',colnames(cors_all_pd) ) ]
cors_all_pd[,'tau_asyn']
progression_markers <- c('nfl_serum', 'lowput_ratio','tau_ab', 'tau_asyn', 'abeta', 'mean_striatum', 'sft', 'td_pigd', 'HVLTRDLY' )
biochemical_markers_conf<-c( 'abeta','tau_ab', 'tau_asyn', 'ptau_ab', 'lowput_ratio', 'ptau', 'ab_asyn'  , 'hemo', 'asyn', 'CSFSAA') # add these biochemical markers 


measured_cells
#sm$updrs3_score_on
clinical_scales_conf<-c('NP2PTOT', 'updrs3_score', 'moca', 'scopa', 'sft', 'sft_V12', 'RBD_TOT', 'NP3TOT_LOG', 'NP3TOT','updrs3_score_on')
clinical_scales<-c(imaging_variables_diff, scale_vars_diff)
MOFAobject@samples_metadata$Plate<-as.factor(MOFAobject@samples_metadata$Plate)
vars_to_plot=c(clinical_scales,progression_markers ); sel_factors<-get_factors_for_scales(clinical_scales)
all_diff_variables_prog<-c(vars_to_plot, 'AGE', 'SEX', 'PDSTATE', 'PD_MED_USE', 'PDMEDYN', 'SITE', 'Plate',  'NP2PTOT_LOG', 'Usable_Bases_SCALE', 
                       'Neutrophils....', 'Lymphocytes....', 'Neutrophils.Lymphocyte', 
                         'sft_V12', c(colnames(estimations),measured_cells), 
                         'Multimapped....', 'Uniquely.mapped....')
all_diff_variables_prog_conf<-c(biochemical_markers_conf, clinical_scales_conf, 'AGE', 'SEX', 'Plate', 'NP2PTOT_LOG','NP3TOT_LOG','NP3TOT',
                    'Neutrophil.Lymphocyte', c(colnames(estimations),measured_cells) )

biochemical_markers_conf
sm$updrs3_score_on
variables_conf_only_clinical<-c(biochemical_markers_conf, clinical_scales_conf, 'AGE', 'SEX','LEDD', 'NP2PTOT_LOG', 'NP3TOT_LOG', 
'updrs3_score_on', 'updrs3_score_on_LOG', 'NP3TOT')
   
sel_factors_conf<-get_factors_for_scales(variables_conf_only_clinical)

graphics.off()
colnames(estimations)



fname<-'factors_covariates_cells_PD'
N_FACTORS[-c(6)]
plot_covars_mofa(selected_covars = c(colnames(estimations),'NP2PTOT_LOG','moca', 'NP2PTOT_LOG_V10', measured_cells),fname,plot,factors=c(1:N_FACTORS)[-c(6)],
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectPD_sel, alpha=0.05 )



fname<-'factors_covariates_cells_PD_cor'
plot_covars_mofa(selected_covars = c(colnames(estimations),'NP2PTOT_LOG','moca', measured_cells),fname,plot='r',factors=1:N_FACTORS,
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectPD_sel )

fname<-'factors_covariates_cells'
plot_covars_mofa(selected_covars = c(colnames(estimations),measured_cells),fname,plot,factors=1:20,
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobject )
                    
fname<-'factors_covariates_cells_HC'
plot_covars_mofa(selected_covars = c(colnames(estimations),measured_cells, 'Multimapped....', 'Uniquely.mapped....', 'Usable_bases_SCALE'),fname,plot,factors=1:20,
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectHC , alpha=0.01)



fname<-'factors_covariates_only_nonzero_strict_PD'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot,factors=sel_factors,
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectPD_sel )
fname<-'factors_covariates_only_nonzero_strict_PD_np3'
plot_covars_mofa(selected_covars=c(all_diff_variables_prog,colnames(estimations)),fname,plot,factors = sel_factors,
            labels_col=TRUE, MOFAobject_to_plot=MOFAobjectPD_sel )

plot='log_pval'

#'LEDD'
#' variables_conf_only_clinical
variables_conf_only_clinical
sel_factors_conf<-get_factors_for_scales(c('NP2PTOT_LOG','moca', 'scopa', 'NP3TOT_LOG','NP3TOT', 'updrs3_score_on', 'updrs3_score_on_LOG'))
sel_factors_conf
fname<-'factors_covariates_strict_PD_conference'
plot_covars_mofa(selected_covars=all_diff_variables_prog_conf,fname,plot,
                 factors = sel_factors_conf,labels_col=TRUE, MOFAobject_to_plot=MOFAobjectPD_sel, res=300 )


fname<-'factors_covariates_strict_PD_conference_clinical'
plot_covars_mofa(selected_covars=variables_conf_only_clinical,fname,plot='log_pval',
                 factors = sel_factors_conf,labels_col=TRUE, MOFAobject_to_plot=MOFAobjectPD_sel, res=300 )


fname<-'factors_covariates_strict_PD_conference_clinical_cor'
plot_covars_mofa(selected_covars=variables_conf_only_clinical,fname,plot='r',
                 factors = sel_factors_conf,labels_col=TRUE, MOFAobject_to_plot=MOFAobjectPD_sel, res=300 )


fname<-'factors_covariates_only_nonzero_strict_cor_PD_np3'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot='r',factors = sel_factors,labels_col=TRUE, MOFAobject=MOFAobjectPD_sel )
graphics.off()

mofa_multi_to_use
estim_cors<-colnames(estimations)[colnames(estimations) %in% colnames(cors_all_pd)]


cors_estim<-cors_all_pd[c(23),  estim_cors]
cors_estim
colnames(cors_estim)[cors_estim>0]

write_vars_output(MOFAobject,vars_by_factor , factors=sel_factors_conf)

write_vars_output



########################## NEW LOGIC GET ONLY FACTORS WITH DIFF IN V16 ##############
############# AFTER WE ADDED THIS TO MOFA


all_diff<-all_diff_variables[all_diff_variables %in% colnames(cors_all_pd)]
all_clin_clusts<-colnames(cors_all_pd)[grep('clin_clust',colnames(cors_all_pd) )]

# 2. scales + 3. imaging S
vars_to_plot=all_diff_variables; 
all_diff_variables_prog<-c(vars_to_plot,'AGE', 'SEX', 'PDSTATE', 'PD_MED_USE', all_clin_clusts, 'NP2PTOT_LOG', 'NP3TOT_LOG')
sm=MOFAobjectPD_sel@samples_metadata
sel_factors_diff<-get_factors_for_scales(all_diff_variables_prog)

all_diff_variables_prog<-all_diff_variables_prog[all_diff_variables_prog %in% colnames(MOFAobjectPD_sel@samples_metadata)]

all_diff_variables_prog=unique(all_diff_variables_prog, c(colnames(estimations),measured_cells))
all_diff_variables_prog
#### Factors related to the longterm change in scale #####
graphics.off()
sel_factors_diff
fname<-'factors_covariates_only_nonzero_strict_PD_diff'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot,factors = sel_factors_diff,labels_col=TRUE, MOFAobject=MOFAobjectPD_sel )
fname<-'factors_covariates_only_nonzero_strict_cor_PD_diff'
plot_covars_mofa(selected_covars=all_diff_variables_prog,fname,plot='r',factors = sel_factors_diff,labels_col=TRUE, MOFAobject=MOFAobjectPD_sel )
fname<-'factors_covariates_only_nonzero_strict_cor_PD'
plot_covars_mofa(selected_covars=selected_covars2_progression,fname,plot='r',factors=get_factors_for_scales(selected_covars2_progression),labels_col=TRUE, MOFAobject=MOFAobjectPD_sel )

factors=sel_factors
### ALSO PLOT DIFF VARS for the controls 
all_diff<-all_diff_variables[all_diff_variables %in% colnames(cors)]
all_cors_diff<-cors[, all_diff]
sel_factors_diff<-which(rowSums(all_cors_diff)>0)
all_diff_variables_prog
all_diff_variables_prog<-c(all_diff_variables, 'AGE', 'SEX', 'PDSTATE', 'PD_MED_USE')
cors


## For this analysis 

sel_factors_coh<-get_factors_for_scales(c('COHORT', 'CONCOHORT'), cors_all=cors)
sel_factors_coh

#### Factors related to the longterm change in scale #####
fname<-'factors_covariates_only_nonzero_strict_diff'
all_diff_variables_prog
sel_factors_diff
plot_covars_mofa(selected_covars=unique(all_diff_variables_prog),fname,plot,factors = 1:25,labels_col=TRUE, MOFAobject=MOFAobject_sel )
fname<-'factors_covariates_only_nonzero_strict_cor'

samples_metadata(MOFAobject_sel)$COHORT<-as.numeric(samples_metadata(MOFAobject_sel)$COHORT)
plot_covars_mofa(selected_covars=c(colnames(estimations), 'COHORT'),fname,plot='r',factors = sel_factors_coh,
labels_col=TRUE, MOFAobject=MOFAobject_sel )

# PD control one with cells and one with clin scores
#  height = 300+100*length(sel_factors_coh),
fname<-'factors_covariates_only_nonzero_strict_cells'
plot_covars_mofa(selected_covars=c(colnames(estimations), 'COHORT'),fname,plot,factors=sel_factors_coh,labels_col=FALSE,
 MOFAobject=MOFAobject_sel, res=300 )
 

fname<-'factors_covariates_only_nonzero_strict'
selected_covars2_progression
plot_covars_mofa(selected_covars=c(selected_covars2, 'COHORT'),fname,plot,factors=sel_factors_coh,labels_col=FALSE,
 MOFAobject=MOFAobject_sel, res=300 ) 


all_diff_variables_prog_in_cors<-all_diff_variables_prog[all_diff_variables_prog %in% colnames(cors_pearson_pd)]





# Plot 1: some more non motor that we discovered
selected_covars_broad=unique(selected_covars_broad)
fname<-'factors_covariates_broad_PD'
MOFAobjectPD@samples_metadata$Lymphocytes....
plot_covars_mofa(selected_covars_broad[! selected_covars_broad %in% measured_cells  ],fname,plot,sel_factors_conf,labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD_sel  )

fname<-'factors_covariates_broad_PD_all_fs'
names(selected_covars_broad)
plot_covars_mofa(selected_covars_broad[! selected_covars_broad %in% measured_cells  ],fname,plot,1:N_FACTORS,labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD_sel  )

fname<-'factors_covariates_broad_cor_PD'
plot_covars_mofa(selected_covars_broad[! selected_covars_broad %in% measured_cells  ],fname,plot='r',sel_factors_conf,labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD_sel  )


cur_names<-read.csv2(paste0(data_dir, '/ppmi/output/curated_names.csv'))$Variable
cur_names



MOFAobject_nams<-MOFAobject
hist(MOFAobject@samples_metadata[,'DYSKIRAT'])
selected_covars<-selected_covars2
selected_covars_pearson<-selected_covars2[!grepl('con_putamen', selected_covars2) ];

f_to_plot<-names(sel_factors)
ind_to_update<-colnames(MOFAobject_nams@samples_metadata) %in%selected_covars_pearson
colnames(MOFAobject_nams@samples_metadata)[ind_to_update]
MOFAobject_nams@samples_metadata$scopa

### Corelation not log_pval####

imaging_variables_diff
fname<-'factors_covariates_img_cor'
plot_covars_mofa(selected_covars=imaging_variables_diff,fname,plot='r',factors,labels_col=FALSE, MOFAobject=MOFAobject )




### 2. Write the covariates for each factor to files ####
### filter only the ones that are correlated 

#write.csv(covariate_corelations, paste0(outdir, '/covariate_corelations.csv'))
dir.create(paste0(outdir, '/covariates/'))
write.csv(cors_pearson, paste0(outdir, '/covariates/covariate_corelations_pearson.csv'))
for (fx in 1:N_FACTORS){
  sig<-cors[fx,]>1.5
  c1<-cors[fx,][sig];
   c2<-cors_pearson[fx,][sig]
  sig_names<-colnames(sig)[which(t(sig))]

  c3<-format(cbind(sig_names,c1,c2), digits=2); c3<-c3[order(c3[,1], decreasing = TRUE),]

  write.csv(c3, 
            paste0(outdir, '/covariates/',fx, '.csv'))
}

top_covariates<-c()
for (fx in 1:N_FACTORS){
  print(fx)
  sig<-cors_all_pd[fx,]>1.5
  c1<-round(cors_all_pd[fx,][sig], digits=2)
  c2<-round(cors_pearson_pd[fx,][sig], digits=2)
  sig_names<-colnames(sig)[which(t(sig))]

  c3<-cbind(c1,c2); 
  rownames(c3  )<-sig_names

  c3<-as.matrix(c3[order(c3[,1], decreasing = TRUE),])
   print(head(c3, n=7 )   )
   top_covariates<-unique(c(top_covariates,rownames(head(c3, n=6))))
}
length(top_covariates)



# fname<-'factors_covariates_top_covariates_PD_all_fs'
# remove_confounders = c(measured_cells)
# plot_covars_mofa(top_covariates[! top_covariates %in% measured_cells  ],fname,plot,1:N_FACTORS,labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD_sel  )

# fname<-'factors_covariates_top_covariates_PD_cor_all_fs'

# #top_covariates_curated<-top_covariates[top_covariates %in% cur_names]
# plot_covars_mofa(top_covariates[! top_covariates %in% measured_cells  ],fname,plot='r',1:N_FACTORS,labels_col=TRUE, height=1500, MOFAobject=MOFAobjectPD_sel  )



view='proteomics'; factor=6

vps=length(MOFAobject@dimensions$D)

fps= as.numeric(MOFAobject@dimensions$K)
views<-names(MOFAobject@dimensions$D)



###########################################################
#### Weights ####
##### WRITE ALL weights for each factor in one file 


### Actually get only factors with higher variance in RNA
dir.create(paste0(outdir, '/top_weights/'))

T=0.3
MOFAobject@dimensions$M
vps = MOFAobject@dimensions$M
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
    print(paste(view, factor))
    all_weights<-MOFA2::get_weights(MOFAobject_gs,views = view, factors=factor, 
                            as.data.frame =TRUE)
    
    ### get the top highly weighted variables - absolute value
    top<-all_weights[order(abs(all_weights$value), decreasing = TRUE),]
    if (high_vars_by_factor[factor, view]){
      write.table(top,paste0(outdir, '/top_weights/top_weights_vals',factor,'_', view,'.txt'), sep = '\t')
      
    }
    
    

  }
  }

graphics.off()
  p1<-plot_variance_explained(MOFAobject, plot_total = T)[[2]]
  p1<-p1+theme(axis.text.x=element_text(size=16), 
               axis.title.y=element_text(size=16), 
               axis.text.y=element_text(size=16))
  p1
  ggsave(paste0(outdir, 'variance_explained_total','.png'),plot=p1, 
         width = 3, height=3, dpi=300)
#install.packages('psych')

  
  
#### Get all weights and put it in one file
  #1. Collate in a list  - order significant
  #2. Stack the lists - 
  
colnames(estimations)


  
  
### wHICH VARIABLES correlate with which factors 
pos_cors<-cors>0  # which are sig. corelated with any factors 
n_factors_pos=1
positive_cors<-cors[,colSums(pos_cors)>n_factors_pos]


### THRESHOLD: IMPORTANT
### significance-which variables to print--> log10(pval)=4 - sig is anything above 1.3
x_cor_t=2

i=111
#### Factor plots ####

dir.create(file.path(paste0(outdir,'/factor_plots/2D/')), showWarnings = FALSE, recursive = TRUE)

to_remove_regex=''
positive_cors_to_plot<-positive_cors[,!grepl(tolower(to_remove_regex),tolower(colnames(positive_cors))) ]
positive_cors_to_plot<-positive_cors[,colnames(positive_cors) %in% c(clinical_scales, selected_covars2, selected_covars_broad)]
colnames(positive_cors_to_plot)
grepl(positive_cors,names(positive_cors_to_plot))
#install.packages('forcats')
select_factors_manually=TRUE
sel_pos_factors_manual=c(4,8)
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
        pos_factors=sel_pos_factors_manual
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
             
      }
}
graphics.off()

met<-MOFAobjectPD@samples_metadata
colnames(MOFAobject@samples_metadata)[grep('Neutr',colnames(MOFAobject@samples_metadata))]
plot(met[, 'Neutrophils....'], met$NP2PTOT)
plot(met[, 'Neutrophils....'], met$updrs2_score)


met$neutrophils<-met[, 'Neutrophils....']

ggplot(met, aes_string(x='Neutrophil.Lymphocyte', y='NP2PTOT'))+
  geom_point(aes_string(x='Neutrophil.Lymphocyte', y='NP2PTOT'))+
  geom_smooth()

hist(met$Neutrophils....)
is.numeric(met[!is.na(met$updrs3_score),'updrs3_score'])
lmfit<-lm(NP2PTOT_LOG~log(Neutrophil.Lymphocyte)+AGE_SCALED+SEX, data=met[!is.na(met$NP2PTOT),])

summary(lmfit)

######## Specific correlations #########
library(tidyverse)

  
MOFAobject@samples_metadata$td_pigd

factors_to_plot<-c(3,4)
factor_cors<-paste0(format(x_cors[factors_to_plot], digits=2), collapse=',')
factor_cors
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


# #### Load features 
# KNOWN GENES
# fs_metadata<-read.csv(paste0(data_dir,'/ppmi/ppmi_data/features_metadata_genes.csv'))
# colnames(fs_metadata)<-c('f', 'feature_id', 'known')
# colnames(fs_metadata)
# fs_metadata$view='RNA'
# fs_met_to_merge<-fs_metadata[,c('feature_id', 'view', 'known'  )]

# source(paste0(script_dir, 'ppmi/mofa_my_utils.R'))
# p_ws<-plot_top_weights2(MOFAobject_gs,
#                        view = 'RNA',
#                        factor = 14,
#                        nfeatures = 15,     # Top number of features to highlight
#                        scale = F
# )
# p_ws


#### 


modify_metadata_for_plotting<-function(MOFAobject){
  #'
  #' 

    rownames(MOFAobject@expectations$W$proteomics_csf)<-gsub('proteomics_csf', 'c', rownames(MOFAobject@expectations$W$proteomics_csf))
   # rownames(MOFAobject@expectations$W$proteomics_csf)<-gsub('_c', '', rownames(MOFAobject@expectations$W$proteomics_csf))
    rownames(MOFAobject@expectations$W$proteomics_plasma)<-gsub('proteomics_plasma', '_p', rownames(MOFAobject@expectations$W$proteomics_plasma))
    
    
    # convert groups to numeric
    sm<-as.data.frame(MOFAobject@samples_metadata)
   # sm[,grep('age', colnames(sm))]<-sapply(sm[,grep('age', colnames(sm))], as.numeric)

    MOFAobject@samples_metadata[,grep('age', colnames(sm))]<-sapply(sm[,grep('age', colnames(sm))], as.numeric)

  return(MOFAobject)

}
MOFAobject


      #  Shorten the rownames 
view = 'proteomics_plasma'
names_vec =   rownames(MOFAobject@expectations$W[[view]])
names_vec


    modify_prot_names(names_vec, view, conv_uniprot = TRUE)
   

modify_features_for_plotting<-function(MOFAobject){
  # modify gene names 
    ens_ids_full<- features_names(MOFAobject)$RNA
    ens_ids<-gsub('\\..*', '', ens_ids_full)
    features_names(MOFAobject)$RNA<-get_symbols_vector(ens_ids)

    # map uniprot ids to gene symbols 
   
    

    rownames(MOFAobject@expectations$W$proteomics_csf)<-gsub('_proteomics_csf', '', rownames(MOFAobject@expectations$W$proteomics_csf))
    rownames(MOFAobject@expectations$W$proteomics_plasma)<-gsub('_proteomics_plasma', '', rownames(MOFAobject@expectations$W$proteomics_plasma))

    rownames(MOFAobject@expectations$W$proteomics_t_csf)<-gsub('_proteomics_t_csf', '', rownames(MOFAobject@expectations$W$proteomics_t_csf))
    rownames(MOFAobject@expectations$W$proteomics_t_plasma)<-gsub('_proteomics_t_plasma', '', rownames(MOFAobject@expectations$W$proteomics_t_plasma))

  
    uniprot_ids<-rownames(MOFAobject@expectations$W$proteomics_plasma)
    gene_symbols<-get_symbol_from_uniprot(uniprot_ids)
# rownames(MOFAobject@expectations$W$proteomics_plasma)
    rownames(MOFAobject@expectations$W$proteomics_plasma)<-gene_symbols$SYMBOL
    
    # replace uniprot
    rownames(MOFAobject@expectations$W$proteomics_csf)<-get_symbol_from_uniprot(rownames(MOFAobject@expectations$W$proteomics_csf))$SYMBOL

    return(MOFAobject)

}


MOFAobject_gs<-modify_features_for_plotting(MOFAobject)
rownames(MOFAobject_gs@expectations$W$proteomics_csf)
#MOFAobject_gs<-modify_metadata_for_plotting(MOFAobject_gs)

W <- get_weights(MOFAobject_gs, factors = 1, views = 4, 
                 as.data.frame = TRUE)

MOFAobject@dimensions$K
views=names(MOFAobject@data)
#views=


nFeatures=15

plot_mofa_top_weighted_all_views(MOFAobject_gs, nFeatures )
#### 3. Save heatmaps and top weighted feaqtures ####
#round(vars_by_factor[fact,])
# Variance threshold: >0.1%
round(vars_by_factor[c(3,13,22),], digits=2)

# rename because value is too long in the legend
MOFAobject@samples_metadata$CONCOHORT_DEFINITION[MOFAobject@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'
MOFAobject_gs@samples_metadata$CONCOHORT_DEFINITION[MOFAobject_gs@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'

cors_heatmap=cors_all_pd
colnames(cors_all_pd['Factor1',])[cors_all_pd['Factor1',]>2]
cors_all_pd[3, ]
graphics.off()

exclude_vars_rec<-colnames(sm)[grep('REC_ID',colnames(sm))]
exclude_vars= c('LAST_UPDATE_M4', 'INFODT_M4', 'NTEXAMTM', 'REC_ID_moca', 'REC_ID_st', 'REC_ID',exclude_vars_rec)

MOFAobject_hm=MOFAobjectPD
#MOFAobject_hm<-modify_features_for_plotting(MOFAobject_hm)
rownames(MOFAobject_hm@expectations$W$proteomics_t_plasma)
rownames(MOFAobject_hm@expectations$W$proteomics_csf)


#MOFAobject_hm<-modify_metadata_for_plotting(MOFAobject_hm)
# Modify the mofa object 
ii=3
# Settings 
cor_p_T<-0.1
cor_T<-ifelse(run_mofa_complete,  1.5, 1.3 )
cor_T
ii=1
vps


MOFAobject_hm
#MOFAobject_hm@samples_metadata<-MOFAobject_hm@samples_metadata[,!duplicated(colnames(MOFAobject_hm@samples_metadata))]
#i=2



rownames(MOFAobject_hm@expectations$W$proteomics_csf)<-modify_prot_names(rownames(MOFAobject_hm@expectations$W$proteomics_csf), 
  view='proteomics_csf', conv_uniprot = TRUE)

rownames(MOFAobject_hm@expectations$W$proteomics_plasma)<-modify_prot_names(rownames(MOFAobject_hm@expectations$W$proteomics_plasma), 
  view='proteomics_plasma', conv_uniprot = TRUE)
rownames(MOFAobject@expectations$W$proteomics_plasma)



for (i in seq(4,vps)){
  for (ii in seq(1,fps)){
 #   for (ii in seq(1,1)){
    print(paste('Modality', i, 'factor', ii))

    cluster_rows=TRUE;cluster_cols=TRUE
    
    
    
    ###### Heatmaps 
    
    nfs=min(round(MOFAobject@dimensions$D[i]*0.05), 30)
   
    #jpeg(paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs, '_cr_',cluster_rows, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
    # Plot heatmaps for each factor only for miRNA 
    
    var_captured<-round(vars_by_factor[ii,i], digits=2)
    main_t<-paste0('Factor ', ii, ', Variance = ',var_captured, '%')
    #log10(0.005) = 2.3
    #cor_T<- -log10(0.005); cor_p_T<-0.15
    modality=names(MOFAobject_hm@dimensions$D)[i]
    if (names(MOFAobject_hm@dimensions$D)[i]=='proteomics'){modality=paste(TISSUE, modality )  }
    main_t<-paste0('Factor ', ii, ', Mod ',modality, ', Variance = ',var_captured, '%')
    
    
    
    ns<-dim(MOFAobject_hm@samples_metadata)[1]
 
    rel_cors<-cors_heatmap[ii,][,cors_heatmap[ii,]>cor_T ]
    rel_cors$moca
    
    # sig holds the names only 
    cors_sig=names(rel_cors); cors_sig
    # remove some
   

    names(rel_cors)
    FT=0
    if (length(cors_sig)==0){
      cors_sig=c()
      rel_cors
    } else if (length(cors_sig)>15){
      FT=min(c(length(cors_sig), 20))
      
      # rel_cors_ordered<-rel_cors[order(-rel_cors)][1:7]
      order(-as.matrix(rel_cors))

      # order the columns 
       rel_cors_ordered<-rel_cors[ order(-as.matrix(rel_cors))]
    

      #rel_cors_ordered<-rel_cors[order(-rel_cors)]
      
      cors_sig<-names(rel_cors_ordered)
      cors_sig<-cors_sig[!(cors_sig %in% exclude_vars)]; cors_sig
      cors_sig<-cors_sig[!grepl( 'LAST_UPDATE|INFO_DT|TM|DT|ORIG_ENTRY|DATE|REC_ID|PAG_', cors_sig)]
      
      cors_sig<-cors_sig[1:FT]




    }
    
    # Check the numbers in the specific view
    patients_in_hm<-colnames(mofa_multi_to_use[[i]])[colnames(mofa_multi_to_use[[i]]) %in% MOFAobject_hm@samples_metadata$PATNO_EVENT_ID]
    patients_in_hm
    plot_heatmap_flag=TRUE
    cors_sig<-cors_sig[cors_sig %in% colnames(MOFAobject_hm@samples_metadata)]
    cors_sig
    MOFAobject_hm@samples_metadata[,cors_sig]
   
    MOFAobject_hm@samples_metadata[,cors_sig]<-sapply(MOFAobject_hm@samples_metadata[,cors_sig], as.numeric)
    if (length(cors_sig)>1){
      
      apply(is.na(MOFAobject_hm@samples_metadata[,cors_sig]),2,any )
      sm<-MOFAobject_hm@samples_metadata
      sm_hm<-sm[sm$PATNO_EVENT_ID %in% patients_in_hm,]
      cors_sig_non_na<-names(which(!(colSums(is.na(sm_hm[,cors_sig]))>NROW(sm_hm)-5) ))
      
      MOFAobject_hm@samples_metadata<-MOFAobject_hm@samples_metadata %>% mutate(across(cors_sig, ~replace_na(., mean(., na.rm=TRUE))))
      # try to see if variance 0 is problematic
      cors_sig_non_na<-names(which(sapply(sm_hm[,cors_sig_non_na], var)>1))
      cors_sig_non_na
      #cors_sig_non_na<-cors_sig_non_na[!cors_sig_non_na %in% c('Neutrophil.Score', 'SITE_V12', 'symptom5_V10', 'symptom5_BL')]

    }else{
      cors_sig_non_na=cors_sig 
    }
length(    cors_sig_non_na)

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

  MOFAobject_hm



    tryCatch({
      p<-plot_data_heatmap(MOFAobject_hm, 
                          view = views[i], 
                          factor =  ii,  
                          features = nfs,
                          groups = groups, 
                          denoise = denoise,
                          cluster_rows = cluster_rows, 
                          cluster_cols = cluster_cols,
                          show_rownames = TRUE, show_colnames = TRUE,
                          scale = "row",
                          na_col = "grey90",
                          annotation_samples = cors_sig_non_na,
                          main=main_t)
   
    p
    #ggsave(hname, plot=p,height=nfs/2, width=(ns+as.numeric(length(cors_sig_non_na) )) )
        if (run_mofa_complete){
          width=ifelse( length(cors_sig_non_na)> 0,ns/10+6,ns/10+4)
          
        }else{
          width=ifelse( length(cors_sig_non_na)> 0,ns/80+6,ns/80+4)
          
        }
        
        ggsave(hname, plot=p,height=nfs/5+2, width=width, dpi=250) 
 
    
  }, 

    error = function(e) {
      an.error.occured <<- TRUE
    print(e)})
  
  
}
  # top weights
  # concat all 
  
  
}  







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





























