
#install.packages('psych')


#install.packages("remotes")
#remotes::install_github("bioFAM/MOFAdata")
##### 
### this one depends on mofa application to inherit: 
### 1. MOFAobject, 2. factors, 
# 3. clinical variables
# 4. outdirs 
# 5. csf/plasma/untargeted flags

#source('enrichment.R')


library(ggplot2)
#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)

graphics.off()
MOFAobject@features_metadata
features_names(MOFAobject)
#getGene(id = rownames(breast_data$mRNA) , type='hgnc_symbol', mart=ensembl )



print(outdir)

jpeg(paste0(outdir, 'factor_cor','.jpeg'))
plot_factor_cor(MOFAobject)
dev.off()
group=1

### Change mofa names to gene symbols 
### create a new object with rna names 
MOFAobject_gs<-MOFAobject

ens_ids_full<- features_names(MOFAobject)$RNA
ens_ids<-   ens_ids<-gsub('\\..*', '', ens_ids_full)


geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ens_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
length(ens_ids)
geneIDs1
length(geneIDs1$SYMBOL)

new_ids<-geneIDs1[match(ens_ids,geneIDs1$GENEID ),]
## sOME SYMBOLS DOI NOT EXIST SO we only replace the ones that do
not_na_ind<-!is.na(new_ids$SYMBOL)
ens_ids[not_na_ind]<-new_ids$SYMBOL[not_na_ind]

features_names(MOFAobject_gs)$RNA<-ens_ids


MOFAobject_gs@samples_metadata$COHORT_DEFINITION

vars_by_factor_all<-calculate_variance_explained(MOFAobject)
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]
write.table(format(vars_by_factor,digits = 2)
            ,paste0(outdir,'variance_explained.txt'), quote=FALSE)


vars_by_factor>0.1

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained','.png'), width = 4, height=4, dpi=100)

MOFAobject@samples_metadata$SEX

colnames(MOFAobject@samples_metadata)[20:35]


stats<-apply(MOFAobject@samples_metadata, 2,table )
sapply(stats,var)>0

non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )

NROW(non_na_vars)
#### Covariance of factors with metadata 


#correlate_factors_with_covariates(MOFAobject,
#                                  covariates = c("NP1ANXS", "NP1DPRS", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
#                                  plot = "log_pval"
#                                 
#                                  
#)
#dev.off()



## plot only factors that correlate? 
cors<-correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars), 
                                  plot = "log_pval", 
                                  return_data = TRUE
                                  
)
ids_to_plot<-which(apply(cors, 2, sum)>0)



jpeg(paste0(outdir, 'factors_covariates_all','.jpeg'), width = 2000, height=700, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars), 
                                  plot = "log_pval"
                                  
)
dev.off()
graphics.off()

jpeg(paste0(outdir, 'factors_covariates_only_nonzero','.jpeg'), width = 2000, height=700, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars)[ids_to_plot], 
                                  plot = "log_pval"
                                  
)
dev.off()

### filter only the ones that are correlated 

covariate_corelations<-correlate_factors_with_covariates(MOFAobject,
                                                         covariates = colnames(MOFAobject@samples_metadata)[c(6:12,45:70, 90:124)], 
                                                         plot = "log_pval",
                                                         return_data = TRUE
)

covariate_corelations<-correlate_factors_with_covariates(MOFAobject,
                                  covariates = colnames(MOFAobject@samples_metadata)[c(6:12,45:70, 90:124)], 
                                  plot = "log_pval",
                                  return_data = TRUE
)
write.csv(covariate_corelations, paste0(outdir, '/covariate_corelations.csv'))


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
  
  all_weights1<-MOFA2::get_weights(MOFAobject,
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

  
  plot_variance_explained(MOFAobject, plot_total = T)[[2]]
  ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)
#install.packages('psych')

  
  
#### Get all weights and put it in one file
  #1. Collate in a list  - order significant
  #2. Stack the lists - 
  

  
  
  
  
  
### wHICH VARIABLES correlate with which factors 
pos_cors<-cors>0  # which have more than two factors positive 
positive_cors<-cors[,colSums(pos_cors)>1]

for (i in 1:dim(positive_cors)[2]){
  

  names<-colnames(positive_cors)
  x_cors<-positive_cors[,i]
  pos_factors<-names(which(x_cors>0))
  
  # Order by 
  pos_factors<-pos_factors[order(x_cors[pos_factors], decreasing = TRUE)]
  print(paste(i, pos_factors))
  
  
  #TODO: print also other factors combinations
  #combn(pos_factors,2)
  
  fs<-c(pos_factors[1],pos_factors[2])
  color_by<-names[i]
  print(color_by)
   
  plot_factors(MOFAobject, 
               factors = fs, 
               shape_by=color_by,
               color_by = 'group',
                 show_missing = FALSE
  )
  
  fss<-paste(fs,sep='_',collapse='-')
  dir.create(file.path(paste0(outdir,'/factor_plots/')), showWarnings = FALSE)
  
  FNAME<-paste0(outdir,'/factor_plots/', 'plot_factors_variate_2D',fss,'_',color_by,'.png')
  
  
  ggsave(FNAME, width = 4, height=4, dpi=100)
  
  shape_by='NHY'
  #shape_by='AGE_AT_VISIT'
  
  plot_factors(MOFAobject, 
               factors = fs, 
               color_by=color_by,
                shape_by= shape_by,
               show_missing = FALSE
  )
  FNAME<-paste0(outdir,'/factor_plots/group/', 'plot_factors_variate_2D',fss,'_',color_by,'_',shape_by,'.png')
    ggsave(FNAME, width = 4, height=4, dpi=100)
  
}


MOFAobject@samples_metadata$CONCOHORT_DEFINITION

# here find the view for which the variability of the factor maximum
plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = fs[1],  
                  features = 4,
                  sign = "positive",
                  color_by = color_by
) + labs(y="Drug response (cell viability)")



## todo extract row and column for which positive cors are TRUE 
#positive_cors_1<-pos_cors>0
#which(positive_cors_1)

#for (i in 1:length(positive_cors_1)){
#  
#}

plot_factors(MOFAobject, 
             factors = fs, 
             show_missing = FALSE
)


plot_factors(MOFAobject, 
             factors = c('Factor1', 'Factor2'), 
             color_by = color_by,
             show_missing = FALSE
)

# Plot top variables with top factors?  
plot_factors(MOFAobject, 
             factors = c(3,4), 
             dot_size = 2.5, 
             color_by = 'NHY'
             
)
dev.off()

##### Plot molecular signatures in the input data



plot_weights(MOFAobject,
             view = "miRNA",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)






### Age, gender , stage does not discriminate factors
#Conclusion here:# factor 1 correlates with grade

  
for (ii in seq(1,fps)){
  ### Plot factors against a clinical variable 
  cors_sig<-names(which(positive_cors[ii,]>2))
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
        
        
        FNAME<-paste0(outdir,'/factor_plots/', 'plot_factors_variate_1D_',ii, '_',color_by,'.png')
        
        
        ggsave(FNAME, width = 4, height=4, dpi=100)
        
        
    }
  }
}



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



##### plot weights 
library(grid)
library(gridExtra)
v_set=c()
v_set=c()

view='miRNA'
factor=8
plot_top_weights(MOFAobject,
                 view = view,
                 factor = factor,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave(paste0(outdir, 'top_weights_',factor, view,'_','.png'), width =3 , height=4, dpi=100)
dir.create(paste0(outdir, 'top_weights/'))

fps=8
vps
fps
seq(1,fps)
seq(1,vps)



graphics.off()

dir.create(paste0(outdir, '/heatmap/'))

for (i in seq(1,vps)){
  for (ii in seq(1,fps)){
   
    
    ### oNLY SAVE THE ones with high variance
    if (high_vars_by_factor[ii, i]){
      print(c(i,ii))
      
          plot_top_weights(MOFAobject_gs,
                           view = views[i],
                           factor = ii,
                           nfeatures = 10,     # Top number of features to highlight
                           scale = T           # Scale weights from -1 to 1
          )
          
          plot_weights(MOFAobject_gs, 
                       view = views[i], 
                       factor = ii, 
                       nfeatures = 30
          )
         
            
          ggsave(paste0(outdir, 'top_weights/all_weights_','f_', ii,'_',vps[i],'.png'), width = 4, height=4, dpi=100)
      
          
          cluster_rows=TRUE;cluster_cols=TRUE
          
          ###### Heatmaps 
          nfs=20
          print('heatmap')
          #jpeg(paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs, '_cr_',cluster_rows, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
          # Plot heatmaps for each factor only for miRNA 
          
          var_captured<-round(vars_by_factor[ii,i], digits=2)
          main_t<-paste0('Factor ', ii, ', Variance = ',var_captured, '%')
         #p<-plot_data_heatmap(MOFAobject, 
         #                     view = views[i], 
         #                     factor =  ii,  
         #                     features = nfs,
         #                     groups=1,
         #                     denoise = TRUE,
         #                     cluster_rows = cluster_rows, cluster_cols = cluster_cols,
         #                     show_rownames = TRUE, show_colnames = TRUE,
         #                     scale = "row",
         #                     fontsize_number = 5, 
         #                     annotation_samples='NHY',
         #                     main=main_t
         #                     
         #                     
         #                     
         #)
         ##main=textGrob(main_t,gp = gpar(fontsize = 21, fontface = "bold") )
         ##grid.arrange(grobs = list(main, p[[4]]), heights = c(0.05, 1))
         #
         #dev.off()
      #          
          ns<-dim(MOFAobject@samples_metadata)[1]
          cor_T<-4
          rel_cors<-cors[ii,][cors[ii,]>cor_T ]
          rel_cors
          
          cors_sig=names(which(cors[ii,]>cor_T))
          FT=0
          if (length(cors_sig)==0){
            cors_sig=c()
            
          } else if (length(cors_sig)>7){
            FT=5
            rel_cors_ordered<-rel_cors[order(-rel_cors)][1:7]
            cors_sig<-names(rel_cors_ordered)
          }
         # exclude_vars= c('LAST_UPDATE_M4')
          #cors_sig<-cors_sig[!(cors_sig %in% exclude_vars)]
          
          res=100
         jpeg(paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, res, '_cor_', cor_T, 'FT_', FT, '.jpeg'),
              height=60*nfs, width=4*ns, res=200)
               
         p<-plot_data_heatmap(MOFAobject_gs, 
                                    view = views[i], 
                                    factor =  ii,  
                                    features = nfs,
                                    denoise = TRUE,
                                    cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                                    show_rownames = TRUE, show_colnames = TRUE,
                                    scale = "row",
                                    annotation_samples=cors_sig,
                                    main=main_t
                                    
                                    
            )
        
      dev.off()
          }
    # top weights
    # concat all 
    
    
    
  }
  
  
}

MOFAobject@samples_metadata$NHY



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
grid.arrange(arrangeGrob(grobs=list(p1, p2), nrow = 1, top="Main Title"))
do.call('grid.arrange', c(list(p1,p2)) )

dev.off()


plot_data_heatmap(MOFAobject, 
                  view = "miRNA",
                  factor = 2,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 3,  
                  features = 100,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


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



### 
## GENE SET ENRICHMENT! 
## AND reactome gsa enrichment!!


###### GSEA 

#BiocManager::install('AnnotationHub')

#source('enrichment.R')
  
#library(AnnotationHub)


library('MOFAdata')
utils::data(reactomeGS)

head((reactomeGS))



subcategory<- 'CP:KEGG'
subcategory<- 'CP:KEGG'
subcategory<- 'GO:MF'
subcategory<- 'GO:BP'
dir.create(paste0(outdir, '/enrichment/'))
for (subcategory in c('GO:BP' ,'CP:KEGG')){
  
        
        gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
        
        gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))
        rownames(gs)
        
        features_names(MOFAobject)$RNA
        features_names(MOFAobject)$RNA<-sapply(features_names(MOFAobject)$RNA, 
               function(x) {stringr::str_remove(x, '\\..*')}
        )
        
        
        # GSEA on positive weights, with default options
        res.positive <- run_enrichment(MOFAobject, 
                                       feature.sets = reactomeGS, 
                                       view = "RNA",
                                       sign = "positive"
        )
        
        
        
        # GSEA on negative weights, with default options
        res.negative <- run_enrichment(MOFAobject, 
                                       feature.sets = gs, 
                                       view = "RNA",
                                       sign = "negative"
        )
        
        
        
        res.positive <- run_enrichment(MOFAobject, 
                                       feature.sets = gs, 
                                       view = "RNA",
                                       sign = "positive"
        )
        
        
        
          
          
          
        # change to negative and positive
        T=0.05
        extract_order_significant<-function(x, T) {
          # extracy most significant and order 
          
          sign<-x[x<T]
          sign2<-sign[order(sign)]
          print(sign2)
        }
        

        stack_list<-function(i,enrichment_list) {
          
          # Take a list of dataframes and stack them 
          # Add a column called factor which extracts the list counter 
          if (length(enrichment_list)){
            
         
                x=enrichment_list[[i]]
                if (length(x)){
                  tmp<-as.data.frame(x)
                  tmp$path<-rownames(tmp)
                  colnames(tmp)<-'pvals'
                  
                  rownames(tmp)<-NULL
                  f<-names(all_fs_enrichment)[[i]] # EXTRACT the counter to assign factor value in new column
                  tmp$factor=f
                  return(tmp)}
                }
          }
        
        
        
        results_enrich<-res.positive$pval.adj
        all_fs_enrichment<-apply(results_enrich, 2 , extract_order_significant, T=T)
        all_fs_unlisted<-sapply(seq(1:length(all_fs_enrichment)), stack_list, enrichment_list=all_fs_enrichment)
        all_fs_merged1<-do.call(rbind, all_fs_unlisted )
        
        write.csv(all_fs_merged1,paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory),  T, '_enrichment_positive_pvals_no_f.csv' ))
        
        all_fs_merged1
        results_enrich<-res.negative$pval.adj
        all_fs_enrichment<-apply(results_enrich, 2 , extract_order_significant,  T=T)
        if (length(all_fs_enrichment)>0){
          all_fs_unlisted<-sapply(seq(1:length(all_fs_enrichment)), stack_list, enrichment_list=all_fs_enrichment)
          all_fs_merged2<-do.call(rbind, all_fs_unlisted )
          
          write.csv(all_fs_merged2,paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory),  T, '_enrichment_negative_pvals_no_f.csv' ))
          #all_fs_merged2
          all_fs_merged2[str_detect(all_fs_merged2[,2], 'PARKINSON'),'factor']
          
        }
       
        ##### which factor is related to parkinsons disease in KEGG
        ### PROBLEM: this is based on RNA only!!! 
    
        
}





# Make enrichment plots for all factors 
# threshold on p value to zoom in 
jpeg(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.positive, 
                        alpha=0.5, 
                        cap=0.0005,
                          colnames=TRUE)
dev.off()

plot_enrichment_heatmap(res.positive$sigPathways, 
                        alpha=0.5, cap=0.0005)

#ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)


jpeg(paste0(outdir,'Enrichment_heatmap_negative','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.negative, 
                        alpha=0.5, cap=0.0005)
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


# Prepare data
# Predict EORTC.risk with factor 1,2 only!
df <- as.data.frame(get_factors(MOFAobject, factors=c(3,4))[[1]])

# Train the model for IGHV
y_predict='NHY'
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
round(importance(model.y), 2)



### Plot predictions
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "NHY",
                  shape_by = "NHY",
                  dot_size = 2.5,
                  show_missing = T
)


cbind(MOFAobject@samples_metadata$sample,MOFAobject@samples_metadata$COHORT_DEFINITION)
MOFAobject@samples_metadata[MOFAobject@samples_metadata$PATNO=='3156',]




