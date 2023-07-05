
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




library(ggplot2)
#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)

graphics.off()

jpeg(paste0(outdir, 'factor_cor','.jpeg'))
plot_factor_cor(MOFAobject)
dev.off()
group=1

### Change mofa names to gene symbols 
### create a new object with rna names 
MOFAobject_gs<-MOFAobject

ens_ids_full<- features_names(MOFAobject)$RNA
ens_ids<-   ens_ids<-gsub('\\..*', '', ens_ids_full)


library(ensembldb)
#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)

## Making a "short cut"
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ens_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new_ids<-geneIDs1[match(ens_ids,geneIDs1$GENEID ),]
not_na_ind<-!is.na(new_ids$SYMBOL)
ens_ids[not_na_ind]<-new_ids$SYMBOL[not_na_ind]
features_names(MOFAobject_gs)$RNA<-ens_ids
MOFAobject_gs@samples_metadata$COHORT_DEFINITION

vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]
write.table(format(vars_by_factor,digits = 2)
            ,paste0(outdir,'variance_explained.txt'), quote=FALSE)

p3<-plot_variance_explained(MOFAobject, max_r2=20)+
  theme(axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))
ggsave(paste0(outdir, 'variance_explained','.png'), plot=p3,
       width = 5, height=N_FACTORS/2,
       dpi=300)

MOFAobject@samples_metadata$SEX

plot_data_overview(MOFAobject)+
  theme(axis.title.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        text  = element_text(size=16))
outdir
ggsave(paste0(outdir, 'data_overview.jpeg'), dpi=300, 
       width=4, height=3)



stats<-apply(MOFAobject@samples_metadata, 2,table )
sapply(stats,var)>0

non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )

NROW(non_na_vars)
#### Covariance of factors with metadata 
source('ppmi/mofa_utils.R')

get_top_cors<-function(MOFAobject){
  cors_both<-get_correlations(MOFAobject, names(non_na_vars))
  cors_top<-cors_both[[1]]
  cors_pearson_top<-cors_both[[2]]
  sel_factors<-abs(cors_pearson_top[,'CONCOHORT'])>0.15
  
  round(cors_top[,'CONCOHORT'][sel_factors], digits=2)
  round(cors_pearson_top[,'CONCOHORT'][sel_factors], digits=2)
}


cors_both<-get_correlations(MOFAobject, names(non_na_vars))
cors_pearson=cors_both[[2]]
cors=cors_both[[1]]
cors_all=cors_both[[1]]

cors_pearson
cors[1,][cors[1,]>1.3]
#max(round(cors_pearson[,'CONCOHORT'][sel_factors], digits=2))



######


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
cors[, '']

MOFAobject@samples_metadata$DATSCAN_CAUDATE_L

#### All associations that are not NA 


jpeg(paste0(outdir, 'factors_covariates_all','.jpeg'), width = 2000, height=700, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = non_na_ids_to_plot, 
                                  plot = "log_pval"
                                  
)
dev.off()
graphics.off()

to_remove_covars<-grepl( 'DATE|REC_ID|UPDATE|ORIG_ENTR|INFO', non_na_ids_to_plot)
non_na_ids_to_plot_cleaned<-non_na_ids_to_plot[!to_remove_covars]

jpeg(paste0(outdir, 'factors_covariates_only_nonzero','.jpeg'), width = length(ids_to_plot)*30, height=2000, res=300)
correlate_factors_with_covariates(MOFAobject,
                                  covariates =non_na_ids_to_plot_cleaned, 
                                  plot = "log_pval"
                                  
)
dev.off()



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
jpeg(paste0(outdir, '/factors_covariates_only_nonzero_strict_pval','.jpeg'),width =700+length(ids_to_plot_strict)*20,
     height = 800, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(ids_to_plot_strict_cleaned), 
                                  plot = "log_pval"
                                  
)
dev.off()



names(non_na_vars)[ids_to_plot]
MOFAobject@samples_metadata$SNCA_rs356181

selected_covars[!(selected_covars %in% names(MOFAobject@samples_metadata))]
#### correlations between factors 
cors[,c('stai_state' )]
format(cors_pearson[,c('stai_trait')], digits=2)

cors[,c('NP1RTOT', 'NP1_TOT','NP2PTOT', 'NP2_TOT','NP3TOT'  ,'NP3_TOT','NP4TOT', 'NP4_TOT' , 'SCAU', 'SCAU_TOT', 'RBD_TOT' )]
format(cors_pearson[,c('NP1RTOT', 'NP1_TOT','NP2PTOT', 'NP2_TOT','NP3TOT'  ,'NP3_TOT','NP4TOT', 'NP4_TOT' , 'SCAU', 'SCAU_TOT','RBD_TOT' )], digits=2)


hist(log2(MOFAobject@samples_metadata$td_pigd_on))
MOFAobject@samples_metadata$MCACLCKH
MOFAobject@samples_metadata$moca
MOFAobject@samples_metadata[,c('con_putamen', 'COHORT', 'updrs3_score_on', 'td_pigd_old_on', 'td_pigd_old', 'cogstate', 'moca', 
                               'ptau', 'ess_cat')]

cors_pearson[,c('NP2PTOT', 'NP2_TOT', 'con_putamen')]
cors
cors[,'updrs3_score_on']

cors[,'NP4_TOT']
cors[,'ess_cat']


selected_covars<-c('COHORT', 'AGE_AT_VISIT', 'SEX', 'NP1TOT', 'NP3TOT', 'NP4TOT', 'SCAU', 'PDSTATE')

selected_covars<-c('COHORT', 'AGE', 'SEX','NP1RTOT', 'NP2PTOT','NP3TOT', 'NP4TOT', 'updrs3_score_on', 
                   'NP1_TOT', 'NP2_TOT','NP3_TOT', 'NP4_TOT',
                   'NHY', 'NP3BRADY',
                   'NP3RIGN', 'SCAU5', 'MCATOT',
                   'MCAVFNUM', 'MCACLCKH', 'cogstate','sft' , 'VLTFRUIT', 'ptau', 'abeta', 'ess_cat', 
                   'HVLTRDLY',
                   'PDSTATE', 'NP3RTCON', 
                  'stai_state', 'stai_trait'  ,'STAIAD26', 'NP1ANXS', 'NP3GAIT', 
                   'SCAU7', 'NP3SPCH', 'NP3RISNG', 'NP2EAT', 
                   'NP3RTARU', 'RBD_TOT', 
                  'con_putamen', 
                 'td_pigd_old_on', 'pd_med_use' )
                   #'DYSKIRAT')
# STAIAD
labels_col=c('Disease status', 'AGE', 'SEX','MDS1','MDS2','MDS3', 'MDS4', 'MDS3_ON',
             'MDS1_log','MDS2_log','MDS3_log', 'MDS4_log',
             'Hoehn & Yahr','MDS3-BRADY','MDS3-RIGN',
             'SC-CONSTIP', 'MCA-cognit-tot'  , 
             'MCA-verb fluency','MCA-visuoconstruct', 'cogstate', 'SEM fluency', 'SEM FL FRUIT',
             'ptau','abeta' ,'EPWORTH sleep cat' ,
             'HOP VERB LEARNING-recall', 
             'PDSTATE', 'MDS3-REST TREMOR', 'STAI_STATE', 'STAI_TRAIT', 'STAI-FEEL RESTED', 'MDSI-anxious', 'MDS3-GAIT', 
             'SC-fec incont', 'MDS3-speech prob', 'MDS3-rising', 'MDS2-eat', 'MDS3-TREMOR', 'RBD_TOT', 
             'PUTAMEN', 
             'TD/PIGD dominant', 'Medication use')
           #  'DYSKIRAT')
# 'STAID:ANXIETY_TOT'

selected_covars<-selected_covars[selected_covars %in% colnames(MOFAobject@samples_metadata)]
labels_col<-labels_col[selected_covars %in% colnames(MOFAobject@samples_metadata)]
graphics.off()
MOFAobject_gs2<-MOFAobject
MOFAobject_gs2@samples_metadata[labels_col]<-MOFAobject_gs2@samples_metadata[selected_covars]
MOFAobject@samples_metadata



jpeg(paste0(outdir, 'factors_covariates_only_nonzero_strict','.jpeg'), width = 800+length(selected_covars)*20, height=1200, res=300)
P2<-correlate_factors_with_covariates(MOFAobject,covariates = selected_covars, plot = "log_pval",
                                  labels_col=labels_col )
dev.off()
selected_covars

jpeg(paste0(outdir, 'factors_covariates_only_nonzero_transpose_strict','.jpeg'),
     height = 800+length(selected_covars)*20, width=1200, res=300)
P2<-correlate_factors_with_covariates(MOFAobject_gs2,covariates = labels_col,
                                      plot = "log_pval", transpose=TRUE
                                       )
dev.off()

plot(P2+coord_flip())


ind_re<-which(non_na_ids_to_plot %in% c('DYSKIRAT'))
## this is the othet
# = non_na_ids_to_plot[-ind_re]

jpeg(paste0(outdir, 'factors_covariates_only_nonzero_cor_logpval','.jpeg'), 
     width = 700+length(selected_covars)*20, height=1100, res=300)
correlate_factors_with_covariates(MOFAobject,covariates = selected_covars, 
                                  plot = "log_pval",
                                 # alpha=0.000000001,
                                  col.lim=c(-0.4, 0.4))
dev.off()

cors_pearson[,c('DYSKIRAT')]

MOFAobject_nams<-MOFAobject
colnames(MOFAobject_nams@samples_metadata)
hist(MOFAobject@samples_metadata[,'DYSKIRAT'])


jpeg(paste0(outdir, 'factors_covariates_only_nonzero_strict_cor','.jpeg'), width = 700+length(selected_covars)*20, height=1100, res=300)
correlate_factors_with_covariates(MOFAobject_nams,covariates =selected_covars,
                                  plot = "r", 
                                  col.lim=c(-0.5, 0.5), 
                                  is.cor=FALSE)
dev.off()



### filter only the ones that are correlated 

#write.csv(covariate_corelations, paste0(outdir, '/covariate_corelations.csv'))
write.csv(cors_pearson, paste0(outdir, '/covariate_corelations_pearson.csv'))

view='proteomics'; factor=6

vps=length(MOFAobject@dimensions$D)
fps= as.numeric(MOFAobject@dimensions$K)
fps
views<-names(MOFAobject@dimensions$D)
views


##### WRITE ALL weights for each factor in one file 


### Actually get only factors with higher variance in RNA
dir.create(paste0(outdir, 'top_weights/'))

T=0.3
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


positive_cors_to_plot<-positive_cors[,!grepl(tolower(to_remove_regex),tolower(colnames(positive_cors))) ]
positive_cors_to_plot
grepl(positive_cors,names(positive_cors_to_plot))

for (i in 1:dim(positive_cors_to_plot)[2]){
  #' fix find 2d plots 
  #' i= index of the clinical variable 

  names<-colnames(positive_cors_to_plot)
  x_cors<-positive_cors_to_plot[,i]
  print(names[i])
  ### does the variable relate to two factors? 
  pos_factors<-names(which(x_cors>0))
  pos_factors<-names(which(x_cors>x_cor_t))
  pos_factors
  # Order by 
  pos_factors<-pos_factors[order(x_cors[pos_factors], decreasing = TRUE)]
  print(paste(i, pos_factors))
  if (length(pos_factors)){
  
          #TODO: print also other factors combinations
          #combn(pos_factors,2)
          if (length(pos_factors)>1){fs=c(pos_factors[1],pos_factors[2])}else{fs=pos_factors[1]}
          fs      
          color_by<-names[i] ## OR change: COLOR BY GROUP? color_by='group'
          print(color_by)
          factor_cors<-paste0(format(x_cors[pos_factors], digits=2), collapse=',')

          pf<-plot_factors(MOFAobject, 
                       factors = fs, 
                       shape_by=color_by,
                       color_by = color_by,
                         show_missing = FALSE 
          )
          
          pf=pf+labs(caption=paste0('log10pval = ',factor_cors))
          fss<-paste(fs,sep='_',collapse='-')
          dir.create(file.path(paste0(outdir,'/factor_plots/2D/')), showWarnings = FALSE)
          
          FNAME<-paste0(outdir,'/factor_plots/2D/', 'plot_factors_variate_2D',fss,'_',color_by, x_cor_t,'.png')
          
          
          ggsave(FNAME,plot=pf, width = 4, height=4, dpi=100)
          
          
          ### you can also add a variable if it is also related to these two factors? 
          shape_by='NHY'
          #shape_by='AGE_AT_VISIT'
          fs
          pf=plot_factors(MOFAobject, 
                       factors = fs, 
                       color_by=color_by,
                        shape_by= shape_by,
                       show_missing = FALSE
          )
          pf=pf+labs(caption=paste0('log10pval = ',factor_cors))
          
          FNAME<-paste0(outdir,'/factor_plots/group/2D/', 'plot_factors_variate_2D',fss,'_',color_by,'_',shape_by, x_cor_t,'.png')
            ggsave(FNAME,plot=pf, width = 4, height=4, dpi=100)
  }
}
dev.off()


MOFAobject@samples_metadata$CONCOHORT_DEFINITION

# here find the view for which the variability of the factor maximum
plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = fs[1],  
                  features = 4,
                  sign = "positive",
                  color_by = color_by
) + labs(y=color_by)


type(MOFAobject@samples_metadata$STAIAD3)

## todo extract row and column for which positive cors are TRUE 
#positive_cors_1<-pos_cors>0
#which(positive_cors_1)

#for (i in 1:length(positive_cors_1)){
#  
#}

color_by='COHORT'
plot_factors(MOFAobject, 
             factors =  c(2,4), 
             color_by = color_by,
             
             show_missing = FALSE
)






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
   
    
    ### oNLY SAVE THE ones with high variance
    if (high_vars_by_factor[ii, i]){
      print(c(i,ii))
      cols <- c( "red", 'red')
          p_ws<-plot_top_weights(MOFAobject_gs,
                           view = i,
                           factor = c(ii),
                           nfeatures = 20,     # Top number of features to highlight
                           scale = F
          )
          graphics.off()
          nFeatures=12
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
          #scale_colour_(values=cols) 
       
            
          ggsave(paste0(outdir, 'top_weights/top_weights_','f_', ii,'_',views[i],'.png'), 
                 plot=p_ws, 
                 width = 2, height=nfeatures/4, dpi=300)
          
          plot_weights(MOFAobject_gs, 
                       view = views[i], 
                       factor = ii, 
                       nfeatures = 30
          )
         
            
          ggsave(paste0(outdir, 'top_weights/all_weights_','f_', ii,'_',views[i],'.png'),
                 width = 4, height=4, dpi=300)
      
          

         
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
    nfs=20
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
    #hname<-paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, res, '_cor_', cor_T, 'FT_', FT, '.jpeg')
    hname<-paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, '_cor_', cor_T, 'FT_', 
                  FT, 'den_', denoise, '.jpeg')

    
    #View(MOFAobject_gs@samples_metadata[cors_sig_non_na])
    p<-plot_data_heatmap(MOFAobject_gs, 
                         view = views[i], 
                         factor =  ii,  
                         features = nfs,
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
        nfs=20
        
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

library('MOFAdata')
utils::data(reactomeGS)

head((reactomeGS))

## TODO: if enrichment is already run then just load results
## load res.positive to be used in the next script

subcategory<- 'CP:KEGG'
subcategory<- 'CP:KEGG'
subcategory<- 'GO:MF'
subcategory<- 'GO:BP'
dir.create(paste0(outdir, '/enrichment/'))
#for (subcategory in c('GO:BP' ,'CP:KEGG')){

mode='proteomics'
mode='RNA'

  for (subcategory in c('GO:BP', 'GO:MF' )){
        if (mode=='proteomics'){
          gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')
          
        }else{
          gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
          
        }
        
        gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))
        colnames(gs)
        
        
        features_names(MOFAobject)$RNA
        features_names(MOFAobject)$RNA<-sapply(features_names(MOFAobject)$RNA, 
               function(x) {stringr::str_remove(x, '\\..*')}
        )
        # GSEA on negative weights, with default options
        res.negative <- run_enrichment(MOFAobject, 
                                       feature.sets = gs, 
                                       view = mode,
                                       sign = "negative"
        )
        
        res.positive <- run_enrichment(MOFAobject, 
                                       feature.sets = gs, 
                                       view = mode,
                                       sign = "positive"
        )
        
        
        
        
        ## TODO: create a function to do for both positive and negative 
        #
        write_enrich<-function(res, sign_mode){
          results_enrich<-res$pval.adj
          all_fs_merged2<-reshape::melt(results_enrich)
          #all_fs_merged2<-all_fs_merged2[all_fs_merged2$value<T,]
          all_fs_merged2<-all_fs_merged2[with(all_fs_merged2, order(X2, value)),]# order 
          
          neg_file<-paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                           mode, '_enrichment', sign_mode)
          write.csv(format(all_fs_merged2, digits=3),paste0(neg_file,  '.csv' ))
          T=0.05
          all_fs_merged2=all_fs_merged2[ all_fs_merged2$value<T,]
          write.csv(format(all_fs_merged2, digits=3),paste0(neg_file, '_', T,  '.csv' ))
          saveRDS(res.negative,paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', sign_mode ))
          
        }
        
        write_enrich(res.negative, sign_mode='negative')
        write_enrich(res.positive, sign_mode='positive')
        
        
        
        
       
        ##### which factor is related to parkinsons disease in KEGG
        ### PROBLEM: this is based on RNA only!!! 
    
        
}





# Make enrichment plots for all factors 
# threshold on p value to zoom in 
jpeg(paste0(outdir,'/enrichment/Enrichment_heatmap_positive','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.positive, 
                        alpha=0.5, 
                        cap=0.0005,
                          colnames=TRUE)
dev.off()

plot_enrichment_heatmap(res.positive$sigPathways, 
                        alpha=0.5, cap=0.0005)

#ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)


jpeg(paste0(outdir,'/enrichment/Enrichment_heatmap_negative','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.negative, 
                        alpha=0.5, cap=0.00000000005 
                          )


dev.off()
#ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


F3<-res.positive$pval.adj[,'Factor3']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:20]

F3<-res.negative$pval.adj[,'Factor6']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:10]



# Positive Factor 1: Hypoxia, oxygen dependent etc.
names(which(res.positive$pval.adj[,'Factor1']<0.0005))[1:10]
names(which(res.negative$pval.adj[,'Factor1']<0.05))

# Positive factor 2: citric acid cycle, meabolism, L13a-mediated, respratory
names(which(res.positive$pval.adj[,'Factor2']<0.00000005))[1:10]
names(which(res.negative$pval.adj[,'Factor2']<0.05))


# Positive factor 3: chromatin modifying enzymes, regulation of tp53
names(which(res.positive$pval.adj[,'Factor3']<0.05))[1:10]
names(which(res.negative$pval.adj[,'Factor3']<0.05))

#  
names(which(res.positive$pval.adj[,'Factor4']<0.05))[1:5]
names(which(res.negative$pval.adj[,'Factor4']<0.05))[1:5]






##### Make predictions 
## Here we select the top factors associated with the clinical variable of interest. 
## then we can choose the top variables of those factors 

## Prediction of clinical subgroups 
## Predict the EORTC.risk

suppressPackageStartupMessages(library(randomForest))

# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(4))[[1]])
df


# Train the model for eortc.risk
df$EORTC.risk <- as.factor(MOFAobject@samples_metadata$EORTC.risk)
model.EORTC.risk <- randomForest(EORTC.risk ~ ., data=df, ntree=10)

# Do predictions
MOFAobject@samples_metadata$EORTC.risk.pred <- stats::predict(model.EORTC.risk, df)
MOFAobject@samples_metadata$EORTC.risk.pred

# Assess performance 
predicted<-as.factor(MOFAobject@samples_metadata$EORTC.risk.pred)
actual = MOFAobject@samples_metadata$EORTC.risk
confusion_mat = as.matrix(table(actual, predicted )) 
print(confusion_mat)

## Show "importance" of variables: higher value mean more important:
round(importance(model.EORTC.risk), 2)
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




