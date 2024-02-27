


# WARNING DO NOT CHANGE THIS, MOVE TO A FUNCTION!! ## 
### FIRST LOAD required files  ####


source(paste0(script_dir, '/ppmi/utils.R'))
knitr_mode<-isTRUE(getOption('knitr.in.progress'))

load_all_se<-function(){
  #'
  #' @param
  #' @return se_rnas, se_mirnas 


      process_mirnas = TRUE; # reload mirs !!  # DO NOT CHANGE THIS!! 


      source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
      if (!base::exists(quote(se_mirs))){
        se_mirs=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
      # se_mirs_norm=load_se_all_visits(input_file = input_file_mirs, combined=combined_bl_log);

      }


      # hist(head(assay(se_mirs_norm),10)[10:50])

      process_mirnas=FALSE
      source(paste0(script_dir, '/ppmi/config.R'))
      #source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))
      if (!base::exists(quote(se_rnas))){

        se_rnas=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
      }

      return(list(se_rnas, se_mirs))

}


all_se_list<-load_all_se()
se_rnas = all_se_list[[1]]
se_mirs = all_se_list[[2]]

## 1. get Summarized Experiment with metrics from all time points 
## 2. Run deseq 
## 3. enrichment 

MOFAobject_clusts=MOFAobject_sel # take it from the clustering of the last visit only 


## SETUP the script parameters here ####
# 1. select visit, 2. process mirs 
# TODO: make function to load for rnas and mirnas separately
# edit this one 

# do not remove these because in the loading..
#process_mirnas= FALSE



if (process_mirnas){
  se_sel = se_mirs
  prefix='mirnas_'
}else{
  se_sel = se_rnas
  prefix='rnas_'
}

view=ifelse(process_mirnas, 'miRNA', 'RNA');view



### Decide on the parameters settings 
# Set the outdirectory 

#y_clust='scopa'
# DIFF_VAR='NP2PTOT_LOG'
#DIFF_VAR='moca'

y_clust = DIFF_VAR

clust_name=paste0(y_clust, '_clust')




# MOFAobject_clusts<-MOFAobjectPD
# Obtain clustering from mofa
se_clusters<-filter_se(se_sel, VISIT=VISIT_COMP, sel_coh = sel_coh, sel_sub_coh = sel_ps)
se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust );
se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)






# MOFAobject_clusts<-MOFAobjectPD
# Obtain clustering from mofa
se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust );
se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)




## Outputs 
# 1. DE files 
# 2. Venns 
# 3. Volcano plot
# 4. Enrichment analysis 


deseq_all_groups <- vector("list", length = 3);
se_filt_all<- vector("list", length = 3);

# Correct for blood cell proportions of neutrophils and lymphocytes 

nclusts = length(table(se_clusters$kmeans_grouping));nclusts
#cluster_params_dir<-paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/')

fact<-get_factors_for_metric(DIFF_VAR); fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors)) 
cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );
cluster_params_dir

cd <- colData(se_clusters)
colData(se_clusters)[cd$INEXPAGE%in%'INEXHC','kmeans_grouping']<-'HC'



se_clusters$kmeans_grouping<-as.factor(se_clusters$kmeans_grouping)




# TODO: include corrected for cell counts  and uncorrected formula
#if (process_mirnas){
  # TODO: try site? and lymphocytes too? 
#  formula_deseq = '~AGE_SCALED+SEX+kmeans_grouping'

# apply this both for mirs and rnas 

variates_to_p<-get_covariates_cells(y_clust, thresh=0)
variates_to_correct<-variates_to_p[!variates_to_p %in% c('Lymphocytes....', 'Neutrophils....', 'Neutrophil.Score')]
variates_to_correct_s<-paste(variates_to_correct, collapse='+')


# add short? 
# variates_to_correct_string



if(!any(colnames(estimations)%in% colnames(colData(se_clusters)))){
  # if cols of estimated cells not inside add them
  # should probably add in the metadata
  #estimations[match(colData(se_clusters)$PATNO_EVENT_ID, rownames(estimations)),]

  estimations_matched<-estimations[match(colData(se_clusters)$PATNO_EVENT_ID, rownames(estimations)),]
  colData(se_clusters)<-cbind(colData(se_clusters),estimations_matched  )

}

# SITE is not required because it is probably explained 
#Error in checkFullRank(modelMatrix) : 
##  the model matrix is not full rank, so the model cannot be fit as specified.
##  One or more variables or interaction terms in the design formula are linear
##  combinations of the others and must be removed


   if (length(variates_to_correct_s)>1){
            variates_to_correct_s<-paste0(variates_to_correct_s, collapse='+')
          }

get_deseq_formula<-function(process_mirnas, cell_corr_deseq,variates_to_correct_s){
      #' get the designformula for deseq based on the type of rnas/mirnas,
      #'  and whether to correct for cell types
      #' @param cell_corr_deseq: correct for cell types
      #' @param formula_deseq_format: correct for neutrophils or all ('n' or 'all')
      formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+COHORT' # basic formula if no cell coreection

          if (process_mirnas){
            # mirnas 
            if (cell_corr_deseq) {

              if ((formula_deseq_format)=='n'){
                formula_deseq = '~AGE_SCALED+SEX+Usable_Bases_SCALE+Neutrophils.LD+COHORT'
              }else{
                formula_deseq = paste0('~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+', variates_to_correct_s,'+COHORT')

              } 
            }
          }else {
            # rnas
            if ((formula_deseq_format)=='n'){
                formula_deseq = '~AGE_SCALED+SEX+Usable_Bases_SCALE+Neutrophils.LD+COHORT'
              }else{
                formula_deseq = paste0('~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+', variates_to_correct_s,'+COHORT')
              }  
           }
  return(formula_deseq)
}


formula_deseq<-get_deseq_formula(process_mirnas, cell_corr_deseq, variates_to_correct_s)



deseq_all <- vector("list", length = 0) # holds the significant gene/mirs ids only for each cluster
deseq_all_names<-vector("list", length = 0)
gene_lists<-vector("list", length = 0)
gse_all<-vector("list", length = 0)

### se_filt_all: list to hold the se
### deseq_all_groups: list to hold the deseq results 
### deseq_significant_all_groups: list to hold significant 
#formula_deseq_format
# if not cell corr ensure that the formula is empty
if (!cell_corr_deseq){
  formula_deseq_format=''
  formula_deseq_format<-''
  print(formula_deseq_format)
}

deseq_params_all<-paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr_deseq))
deseq_params<-paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr_deseq),  '/',VISIT_COMP, '/', formula_deseq_format, '/')
dir.create(paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr_deseq),  '/',VISIT_COMP, '/', formula_deseq_format, '/'))
dir.create(deseq_params, recursive = TRUE)
dir.create(cluster_params_dir, recursive = TRUE)

fname_venn=paste0(deseq_params,'/', prefix , 'min_',min.count,'venn.png');fname_venn



clusters_indices= list( c('1', 'HC'),c('2', 'HC'),c('3', 'HC'), c('1','2','3', 'HC'), c('1','2'), c('1','3'),c('3','2'), c('3', '1') )
# the order for deseq levels is first controls, second disease
deseq_order= list( c('1', '2'),c('1', '2'),c('1', '2') ,c('1', '2'), c('1','2'), c('1','3'),c('3','2'), c('3', '1') )

clusters_names<-c(1,2,3,'1_2_3', '1_2', '1_3', '3_2' , '3_1')


  se_filt_all<-lapply(clusters_indices, function(cluster_id){ 

    se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id)]} )


names(se_filt_all)<-clusters_names

se_filt_all<-lapply(se_filt_all, function(se_clust){
    if (length(unique(na.omit(se_clust$COHORT)))==1){
    se_clust$COHORT<-se_clust$kmeans_grouping}
    return(se_clust)
    
  }
)

# todo remove the ones below
se_filt_all[['1_2']]$COHORT<-se_filt_all[['1_2']]$kmeans_grouping




cluster_id_num = 1
#dir.create(deseq_params, recursive = TRUE)
for (cluster_id_num in 1:length(clusters_names)){

  cluster_id_index=clusters_indices[[cluster_id_num]]
  cluster_id_name=clusters_names[[cluster_id_num]]
 

  de_file<-paste0(deseq_params, '/', prefix, 'de_cluster_', cluster_id_name , '.csv')
  de_params_file<-paste0(deseq_params, '/', prefix, make.names(formula_deseq))

   if (file.exists(de_file)){
  #if (FALSE)
    # if de file exists load it - unfiltered de results file
    deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )

  }else{
    print(de_file)
    # else run the deseq with the design formula specified 
        se_clust<-se_filt_all[[cluster_id_name]]
        deseq2ResDF = deseq_by_group(se_clust, formula_deseq, min.count=min.count, contrast_order=deseq_order[[cluster_id_num]])
        #deseq2ResDF$log2FoldChange


        # pipeline
        # 1. run deseq
        # 2. save results 
        # 3. plot Volcano 
        # 4. run enrichment 

        deseq_all_groups[[cluster_id_name]]<-deseq2ResDF
        if (!process_mirnas){
          # get symbols for RNA only 
          deseq2ResDF$GENE_SYMBOL<-get_symbols_vector(gsub('\\..*', '',rownames(deseq2ResDF))) 
        }
        write.csv(deseq2ResDF, de_file, row.names=TRUE)
        file.create(de_params_file)


        # Volcano
        pvol<-plotVolcano(  deseq2ResDF, se_clust, title=paste0('Cluster ', cluster_id_name), xlim=c(-1.1,1.1),
            lab=deseq2ResDF$GENE_SYMBOL)

      fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT_COMP,'.jpeg')
      fname<-paste0(deseq_params, '/Volcano_', prefix, VISIT_COMP,'_cluster_', cluster_id_name, '.jpeg')
      ggsave(fname,pvol, width=9,height=12, dpi=300)



        # Also write the parameters of deseq

  }

 

  deseq_all_groups[[cluster_id_name]]<-deseq2ResDF
  deseq_all[[cluster_id_name]]<-deseq2ResDF[deseq2ResDF$mofa_sign %in% 'Significant',] # holds the significant only
} 
# Save and load # Rrename ens id.*

de_file
deseq_all_names <- lapply(deseq_all, function(x){return(  gsub('\\..*', '',rownames(x))   )  })
names(deseq_all_names) <-names(deseq_all)

deseq_all_names


#### 1. Venn from significant 
length(deseq_all_names)
head(deseq_all_names,3)
create_venn(venn_list = head(deseq_all_names,3), fname_venn =fname_venn,
            main =paste0( ' DE molecules for each molecular cluster' ) )

graphics.off()
 # TODO: 

########### 



#### 2. Venn from significant in top of factor

## 1. get top of factor / or all higly variable genes input into MOFA
## 2. intersect with DE 
sel_factor=4
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=view, factors = sel_factor, top_fr = 0.1))
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=view, factors = sel_factor, top_fr = 1))
length(top10)
# intersect with the top factors 
deseq_all_top<-lapply(deseq_all_names, function(x) intersect(x,top10 ) )
deseq_all_top


names(deseq_all_top)
length(deseq_params)
#fname_venn=paste0(deseq_params, prefix, 'venn_de_per_group_deseq', 'top_f', sel_factor ,  '.png')
#create_venn(venn_list = deseq_all_top, fname_venn =fname_venn,main =paste0( ' DE molecules for each molecular cluster AND highly variable' ))





order_by_metric='log2pval';order_by_metric_s='log2p'
order_by_metric='log2FoldChange'; order_by_metric_s='log2FC'


pvalueCutoff_sig=0.05
enrich_params<-paste0(ONT, order_by_metric_s)
dir.create(paste0(deseq_params, '/enr/'))

deseq_params
enrich_compare_path=paste0(deseq_params, '/enr/', prefix, enrich_params, 'comp')

# run enrichment only for RNAs 

if (!process_mirnas){

gse_all_clusters=list()

force_gse=FALSE # SET to true to rerun plots 
#clusters_names
 for (cluster_id_num in 1:length(clusters_names)){

      print(paste('cluster:',cluster_id_num))
      cluster_id_index=clusters_indices[[cluster_id_num]]
      cluster_id_name=clusters_names[[cluster_id_num]]
       
       # parse the deseq results  for the analysis
      deseq2ResDF = deseq_all_groups[[cluster_id_name]]
     # parse the gene list from the deseq results 
      gene_list1<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )
      names(gene_list1)<-gsub('\\..*', '',names(gene_list1))

      # load the gene lists separately for cluster compare  
      # TODO: move to another loop to save this intermediate output 
      gene_lists[[cluster_id_name]]<-gene_list1 

      # get the gene list and run ! - maybe it is not needed if we parse the cluster compare....
      results_file_cluster=paste0(deseq_params, '/enr/', prefix, enrich_params, 'cl', cluster_id_name)
      gse_file<-paste0(results_file_cluster, '.Rds')

     # if (TRUE){
      if (!file.exists(gse_file) | force_gse){

          gse1<-run_enrich_per_cluster(deseq2ResDF, results_file_cluster,N_DOT=20, 
          N_EMAP=30 , N_NET=10, force_gse=FALSE)
          saveRDS(gse1,gse_file )
          gse_all_clusters[[cluster_id_name]]<- gse1

      }else{
        gse_all_clusters[[cluster_id_name]]<- loadRDS(gse_file )
      }
    }







gse_sig_all_clusters<-lapply(gse_all_clusters, function(x) {  x@result$ID[x@result$p.adjust<0.05] } )
fname_venn=paste0(deseq_params, '/',prefix, 'venn_de_per_group_enrich.png')

#if (length(gse_sig_all_clusters)>0){
#  create_venn(venn_list = gse_sig_all_clusters, fname_venn =fname_venn,
#main =paste0( ' DE pathways for each molecular cluster' ) )
#}




force_compare=FALSE
# Rerun them all together so that they are in one file for comparisons 
cluster_pairs<-c('1','1_2', '1_3', '2', '3')# Unique in 1: not in 2,3 and also 
cluster_pairs<-c('1','2', '3') 
cluster_pairs<-c('1','1_2', '1_3', '2', '3','3_1', '3_2','1_2_3' ) 


gene_lists_compare<-gene_lists[cluster_pairs]
enrich_compare_path_all<-paste0(deseq_params, '/enr/', prefix, enrich_params, 'comp')
enrich_compare_path <- paste0( enrich_compare_path_all, paste0(names(gene_lists_compare), collapse='-'))
gse_compare_file<-paste0(enrich_compare_path, '.Rds')
gse_compare_file
N_EMAP_COMP=60

if (!file.exists(gse_compare_file) | force_compare){
    print('Running comparisons')
    # Run cluster compare by cluster - it does not need the separate files only the gene lists 
    # for each cluster for ONE time point 
    geneClusters=gene_lists_compare


    gse_compare<-compareCluster(geneClusters = geneClusters , 
                                fun = "gseGO", 
                                OrgDb='org.Hs.eg.db', 
                              ont=ONT, 
                              keyType = 'ENSEMBL', 
                              pvalueCutoff=1) 


    ### RUN SCRIPT compare
    plot_enrich_compare(gse_compare,paste0(enrich_compare_path), N_EMAP = N_EMAP_COMP, N_DOT=8,N_DOT_U = 20, pvalueCutoff_sig = 0.05)
    saveRDS(gse_compare, gse_compare_file)
}else{
    gse_compare<-loadRDS(gse_compare_file)
}


force_compare_all_options=TRUE
if (force_compare_all_options | !knitr_mode){

    
# Subset the comparisons and plot different versions 
filt_clusts_all<-list( c('1','1_2', '1_3', '2', '3'), # unique in 1
                        c('3','3_1', '3_2', '1', '2'), # unique in 3
                        c('1','2','3'), 
                        c('1','2','3', '1_2_3')) # just compare with controls 

lapply(filt_clusts_all, function(filt_clusts){
        # subset and compare 
        # subset 
        enrich_compare_path <- paste0( enrich_compare_path_all , paste0(filt_clusts, collapse='-'))
        pvalueCutoff_sig = 0.05
        #   
        # Check if already run eg. if emap is there

        enrich_compare_path
        gse_compare_sub<-gse_compare %>% 
          dplyr::filter(Cluster %in% filt_clusts)


          plot_enrich_compare(gse_compare_sub,paste0(enrich_compare_path), N_EMAP = 60,
          N_DOT=8,N_DOT_U = 20, pvalueCutoff_sig = pvalueCutoff_sig)
        }
      
)

}



}





































 # nolint




















































































