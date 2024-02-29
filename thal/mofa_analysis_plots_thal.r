


# get correlatins for selected vars 





######### Cluster patients with clinical trajectory factors ####
### get metrics to test how many clusters to use ####
colnames(sm)
sm$Biological.group
#BiocManager::install('M3C')
#library(M3C)
library('cluster')


y='Biological.group'

sm<-samples_metadata(MOFAobject)


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

if (cluster_samples_mofa){
  if (length(sel_coh)>1){
    #for (k_centers_m in c(6)){
    for (k_centers_m_try in c(3,2)){
      clusters <- cluster_samples(MOFAobject, k=k_centers_m_try, factors=c(1:15))
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

MOFAobjectPatients

####
# FACTORS COVARIATES 
met<-samples_metadata(MOFAobjectPatients)
MOFAobjectPD_sel@samples_metadata$Biological.group<-as.factor(MOFAobjectPD_sel@samples_metadata$Biological.group)
graphics.off()
fname<-'factors_covariates_all'
plot_covars_mofa(selected_covars = colnames(sm),fname,plot,factors=c(1:N_FACTORS)[-c(6)],
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobject, alpha=0.05 )

MOFAobjectPatients@samples_metadata$Biological.group<-as.numeric(as.factor(MOFAobjectPatients@samples_metadata$Biological.group))

fname<-'factors_covariates_patients'
plot_covars_mofa(selected_covars = colnames(sm),fname,plot,factors=c(1:N_FACTORS),
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectPatients, alpha=0.05 )


fname<-'factors_covariates_patients_cor'
plot_covars_mofa(selected_covars = colnames(sm),fname,plot='r',factors=c(1:N_FACTORS),
                    labels_col=FALSE, MOFAobject_to_plot =MOFAobjectPatients, alpha=0.05 )





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


get_factors_for_scales<-function(vars_to_plot, cors_all=cors_all_pd){
  #''
  #' @param
  #'
  #'
  vars_to_plot<-vars_to_plot[vars_to_plot %in% colnames(cors_all)]
    all_cors_diff<-cors_all[, vars_to_plot]
  sel_factors<-which(rowSums(all_cors_diff)>0)
  return(sel_factors)
  
}

MOFAobject@samples_metadata$tau_asyn




#3 categorical
sm$Biological.group
clinical_scales_cat<-c('Biological.group', 'SEX')
correlate_factors_categorical('Biological.group')


  cors_kruskal<-sapply(clinical_scales_cat, correlate_factors_categorical, MOFAobject=MOFAobject)
  rownames(cors_kruskal)<-1:N_FACTORS
  kw_cors<-format(cors_kruskal, digits=1)<0.05
  factors_cor_with_clusters<-kw_cors
  #library(pheatmap)
  jpeg(paste0(outdir, 'heatmap_categorial_covariates.jpeg'), units='in', width=4, height=4, res=300)
  pheatmap(-log10(cors_kruskal))
  dev.off()

  cors_kruskal<-sapply(clinical_scales_cat, correlate_factors_categorical, MOFAobject=MOFAobjectPatients)
  rownames(cors_kruskal)<-1:N_FACTORS
  kw_cors<-format(cors_kruskal, digits=1)<0.05
  factors_cor_with_clusters<-kw_cors
  #library(pheatmap)
  jpeg(paste0(outdir, 'heatmap_categorial_covariates_patients.jpeg'), units='in', width=4, height=4, res=300)
  pheatmap(-log10(cors_kruskal))
  dev.off()





variables_clinical<-c(colnames(sm))
colnames(sm)
MOFAobjectPD_sel
sel_factors_conf<-get_factors_for_scales(variables_clinical)

graphics.off()



###########################################################
#### Weights ####
##### WRITE ALL weights for each factor in one file 


### Actually get only factors with higher variance in RNA
dir.create(paste0(outdir, '/top_weights/'))

T=0.3
MOFAobject@dimensions$M
vps = MOFAobject@dimensions$M
MOFAobject@dimensions$M

views = names(MOFAobject@dimensions$D)
views
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
  plot(p1)
  ggsave(paste0(outdir, 'variance_explained_total','.png'),plot=p1, 
         width = 3, height=3, dpi=300)
#install.packages('psych')

  
  
#### Get all weights and put it in one file
  #1. Collate in a list  - order significant
  #2. Stack the lists - 
  
colnames(estimations)


which.min(Z[,5])
  
### wHICH VARIABLES correlate with which factors 
cors<-as.data.frame(cors)
pos_cors<-cors>0  # which are sig. corelated with any factors 
n_factors_pos=1
cors
colSums(pos_cors)
colSums(pos_cors)>1
cors
positive_cors<-cors[,which(colSums(pos_cors)>1)]


### THRESHOLD: IMPORTANT
### significance-which variables to print--> log10(pval)=4 - sig is anything above 1.3
x_cor_t=2




i=111
#### Factor plots ####
#### MANUAL
fs=c(1,4)
total_vols<-colnames(sm)[grep('Total|Number',colnames(sm))]
MOFAobject@samples_metadata[,total_vols] <-sapply(MOFAobject@samples_metadata[,total_vols], as.numeric)
MOFAobjectPatients@samples_metadata$FIB.4<-as.numeric(MOFAobjectPatients@samples_metadata$FIB.4)
MOFAobjectPatients@samples_metadata$APRI<-as.numeric(MOFAobjectPatients@samples_metadata$APRI)
colnames(sm)
MOFAobjectPatients@samples_metadata$APRI
color_by='Biological.group'
fss<-paste(fs,sep='_',collapse='-')
dir.create(file.path(paste0(outdir,'/factor_plots/2D/')), showWarnings = FALSE)
FNAME
FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by,'.png')
     pf<-plot_factors(MOFAobject,  factors = fs, color_by = color_by, show_missing = FALSE )
                      
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)

 FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by, 'only_patients','.png')
pf<-plot_factors(MOFAobjectPatients,  factors = fs, color_by = color_by, show_missing = FALSE )              
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)

 
 FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by,shape_by, 'only_patients','.png')
sm$X2
 color_by = 'FIB.4'; shape_by='Biological.group'
pf<-plot_factors(MOFAobjectPatients,  factors = fs, shape_by=shape_by,color_by = color_by, show_missing = FALSE )              
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)

  color_by = 'APRI'; shape_by='Biological.group'

 FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by,shape_by, 'only_patients','.png')
pf<-plot_factors(MOFAobjectPatients,  factors = fs, shape_by=shape_by,color_by = color_by, show_missing = FALSE )              
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)

   color_by = 'X2'; shape_by='Biological.group'

 FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by,shape_by, 'only_patients','.png')
pf<-plot_factors(MOFAobjectPatients,  factors = fs, shape_by=shape_by,color_by = color_by, show_missing = FALSE )              
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)

  color_by = 'FERITIN'; shape_by='Biological.group'
 FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by,shape_by, 'only_patients','.png')
pf<-plot_factors(MOFAobjectPatients,  factors = fs, shape_by=shape_by,color_by = color_by, show_missing = FALSE )              
 ggsave(FNAME,plot=pf, width = 3, height=3, dpi=200)











to_remove_regex=''
clinical_scales=colnames(sm)
 covars_to_plot<-c(clinical_scales)
positive_cors_to_plot<-positive_cors[,!grepl(tolower(to_remove_regex),tolower(colnames(positive_cors))) ]
positive_cors_to_plot<-positive_cors[,colnames(positive_cors) %in% covars_to_plot]
colnames(positive_cors_to_plot)
grepl(positive_cors,names(positive_cors_to_plot))
#install.packages('forcats')
select_factors_manually=TRUE
sel_pos_factors_manual=c(1:N_FACTORS)[c(4,5)]

sel_pos_factors_manual
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
color_by='Biological.group'

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

tail(MOFAobject_gs@features_metadata$fea)

tail(features_metadata(MOFAobject_gs))


#### 

rownames(MOFAobject_gs@expectations$W$proteomics_csf)<-gsub('proteomics_csf', 'c', rownames(MOFAobject_gs@expectations$W$proteomics_csf))
rownames(MOFAobject_gs@expectations$W$proteomics_csf)<-gsub('_c', '', rownames(MOFAobject_gs@expectations$W$proteomics_csf))
rownames(MOFAobject_gs@expectations$W$proteomics_plasma)<-gsub('proteomics_plasma', '_p', rownames(MOFAobject_gs@expectations$W$proteomics_plasma))



W <- get_weights(MOFAobject_gs, factors = 1, views = 4, 
                 as.data.frame = TRUE)

MOFAobject@dimensions$K
views=names(MOFAobject@data)
#views=

  for (i in seq(1,MOFAobject@dimensions$M)){
  for (ii in seq(1,MOFAobject@dimensions$K)){
   
    nFeatures=20
    
    ### oNLY SAVE THE ones with high variance
   # if (high_vars_by_factor[ii, i]){
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
            #plot_top_weights2
            p_ws<-plot_top_weights(MOFAobject_gs,
                                    view = 'RNA',
                                    factor = ii,
                                    nfeatures = nFeatures,     # Top number of features to highlight
                                    scale = F
            )

          
            print('rna')
            p_ws<-p_ws+ theme(axis.text.y = element_text(face = "italic"))
          }
          width=3
          if (views[i]=='metabolites'){
           p_ws<-p_ws+ theme(axis.text.y = element_text(size=7))
            width=5
           
          }
          
          #scale_colour_(values=cols) 
       
          try(  
          ggsave(paste0(outdir, 'top_weights/top_weights_','f_', ii,'_',views[i],'_',nFeatures,'.png'), 
                 plot=p_ws, 
                 width = width, height=(nFeatures/5)+0.5, dpi=300)
          )
          
          
        #  plot_weights(MOFAobject_gs, 
        #               view = views[i], 
         #              factor = ii, 
        #               nfeatures = 30
         # )
         
            
        #  ggsave(paste0(outdir, 'top_weights/all_weights_','f_', ii,'_',views[i],'_',nFeatures, '.png'),
        #         width = 3, height=nFeatures/4, dpi=300)
      
          

         
          
  }
}




    # top weights
    # concat all 
    
    
 


#### 3. Save heatmaps and top weighted feaqtures ####


# rename because value is too long in the legend
MOFAobject@samples_metadata$CONCOHORT_DEFINITION[MOFAobject@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'
MOFAobject_gs@samples_metadata$CONCOHORT_DEFINITION[MOFAobject_gs@samples_metadata$CONCOHORT==0]<-'non-PD, non-Prod, non-HC'

dim(cors_all_pd)
cors_heatmap=cors_all_pd
cors_all_pd[3, ]
graphics.off()

    exclude_vars= c('LAST_UPDATE_M4', 'INFODT_M4', 'NTEXAMTM', 'REC_ID_moca', 'REC_ID_st')
MOFAobject_hm=MOFAobjectPD

ii=3
fps=N_FACTORS
for (i in seq(1,2)){
  for (ii in seq(1,fps)){
    #print(paste('Modality', i, 'factor', ii))

    cluster_rows=TRUE;cluster_cols=TRUE
    
    
    
    ###### Heatmaps 
    nfs=40
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
    if (run_mofa_complete){
      cor_T<-1.5; cor_p_T<-0.1
      
    }else{
      cor_T<-2; cor_p_T<-0.1
      
    }
    rel_cors<-cors_heatmap[ii,][,cors_heatmap[ii,]>cor_T ]
    rel_cors
    # sig holds the names only 
    cors_sig=names(rel_cors); cors_sig
    FT=0
    if (length(cors_sig)==0){
      cors_sig=c()
      rel_cors
    } else if (length(cors_sig)>15){
      FT=15
      # rel_cors_ordered<-rel_cors[order(-rel_cors)][1:7]
      order(-as.matrix(rel_cors))
       rel_cors_ordered<-rel_cors[ order(-as.matrix(rel_cors)),][1:FT]
      #rel_cors_ordered<-rel_cors[order(-rel_cors)]

      cors_sig<-names(rel_cors_ordered)
    }

    cors_sig<-cors_sig[!(cors_sig %in% exclude_vars)]; cors_sig
    cors_sig<-cors_sig[!grepl( 'LAST_UPDATE|INFO_DT|TM|DT|ORIG_ENTRY|DATE|PAG_', cors_sig)]
    
    plot_heatmap_flag=TRUE
    cors_sig<-cors_sig[cors_sig %in% colnames(MOFAobject_hm@samples_metadata)]
    MOFAobject_hm@samples_metadata[,cors_sig]
    
    #is.na(MOFAobject_gs@samples_metadata[,cors_sig])

    ### if the col contains only NA
    
    #which(cors_sig_non_na=='PDSTATE')
    #cors_sig_non_na=cors_sig_non_na[-3]
    if (length(cors_sig)>1){
      cors_sig_non_na<-names(which( !apply(is.na(MOFAobject_hm@samples_metadata[,cors_sig]),2,any )))
      
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







