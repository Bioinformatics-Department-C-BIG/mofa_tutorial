
DIFF_VAR
# cm_all=list()


### Heatmaps of the log2FC values 
# Input: log2FC from deseq with controls for each gene for each time points 
# Read in deseq results from each timepoint 
# And each cluster 
# need to run after DE tutorial ppmi cluster_compare.R



view = 'RNA'
prefix='rnas_'


view = 'miRNA'
prefix='mirnas_'

#write.csv(results_de, paste0(outdir_s_p, 'results.csv'))
prefix='prot_'
tissue_un_mofa<-ifelse(tissue_un=='Plasma', 'plasma', 'csf')

#' @param tissue is a global name that adjusts for targeted or untargeted 
if (prot_de_mode=='t'){
        view=paste0('proteomics_t_', tolower(TISSUE))
        top_fr=0.025
        tissue=TISSUE
        metric_p<-'adj.P.Val'; T_p=0.05 
         metric_p='P.Value';  T_p=0.001

        metric_p_lfc  = 'logFC'

        sig_only=TRUE
        view_s = paste0(tissue, prot_de_mode_s)



}else if (prot_de_mode=='u'){
        view=paste0('proteomics_', tolower(tissue_un_mofa))
        print(view)
        tissue = tissue_un
        top_fr=0.25
        metric_p<-'adj.P.Val'; T_p=0.05 
        sig_only=FALSE
        view_s = paste0(tissue, prot_de_mode_s)




}

clust_ids<-c('1', '2', '3', '4')
clusters_names_h = c( "1"  ,   "2" ,    "3" ,    "all")
clusters_names =  c( "1"  ,   "2" ,    "3" ,    "all")



if (view %in% c('RNA', 'miRNA')){
       metric_p<-'padj'; 
       metric_p_lfc  = 'log2FoldChange'

       
        clust_ids<-c('1', '2', '3', '1_2_3')
        clusters_names_h = c( "1"  ,   "2" ,    "3" ,    "all")
        top_fr=0.015
        sig_only=TRUE
        view_s = paste0(view)


}
if (view %in% c( 'miRNA')){

        top_fr=0.025
        sig_only=FALSE
        view_s = paste0(view)


}



# TODO: update also for mirnas here 
#cluster_params_dir
# outdir_s_p
fact = get_factors_for_metric(DIFF_VAR)

#fact = fact[fact!=13]
fact
cluster_params_dir<-get_cluster_params_dir(DIFF_VAR)
cluster_params_dir
top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact, view=view, top_fr=top_fr   )
#top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact2, view=view, top_fr=top_fr   )

top_proteins
top_proteins$feature<-gsub(paste0('_',view),'', top_proteins$feature)

top_proteins$feature



cluster_params_dir
get_cluster_de_result_file<-function(outdir_s_p, view, cluster_id){


                    if (grepl('prot', view)){
                        de_file_cluster<-paste0(outdir_s_p, prefix, tissue,'_', prot_de_mode,'_de_cl',cluster_id,  '_results.csv')

                    }else if (view == 'RNA'){
                    # example:    rnas_de_cluster_2
                        de_file_cluster<-paste0(outdir_s_p, prefix, 'de_cluster_',cluster_id,  '.csv')

                    }else if (view == 'miRNA'){
                        de_file_cluster<-paste0(outdir_s_p, prefix, 'de_cluster_',cluster_id,  '.csv')

                    }
                    return(de_file_cluster)
                }




de_sig_all<-list()
VISIT_COMP_SIG = 'V08'


cluster_id = 1
cluster_id = 2
# TODO: this part is in progress- working on making it applicable to RNA modality as well !! 
  for (cluster_id in clust_ids){
        for (VISIT_COMP_SIG in c('BL', 'V06', 'V04', 'V08')){

        
                outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP_SIG, '/' )
                de_prot_file<-paste0(outdir_s_p, prefix, tissue,'_', prot_de_mode,'_de_cl',cluster_id,  '_results.csv')
                de_prot_file = get_cluster_de_result_file(outdir_s_p, view, cluster_id)
                de_results_prot<-read.csv(de_prot_file)
                #print(de_results_prot[,metric_p ])
                
               # view=paste0('proteomics_', tolower(TISSUE))
                # TODO: take the top MOFA proteins from moca 
              
                #}
                colnames(de_results_prot)
                de_results_prot_sig<-de_results_prot[de_results_prot$P.Value<0.01,] # take only the union of all at the end 
                pass_pval<-which(de_results_prot[, metric_p]<0.05)
                de_results_prot_sig<-de_results_prot[pass_pval,] # take only the union of all at the end 
               # print(paste(cluster_id))
              #   print(de_results_prot_sig[de_results_prot_sig$GENE_SYMBOL == 'GZMH',])


                de_sig_all[[cluster_id]]<-de_results_prot_sig$X
                }

  }
de_sig_all_top<-unique(unlist(de_sig_all))
de_sig_all_top

de_sig_all
de_results_prot_sig
de_results_prot_sig$GENE_SYMBOL
# TODO: separate to get top 






get_de_proteins_per_tp<-function(VISIT_COMP, metric_p='logFC', sig_only =FALSE, de_sig_all_top){
        #' 
        #' @param  VISIT_COMP
        #' @param metric_p metric to use to cut 
        #' @param sig_only filter the significant otherwise the top in mofa 

        de_all<-list()
        for (cluster_id in clust_ids){

                outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
                # 
                de_prot_file<-get_cluster_de_result_file(outdir_s_p, view, cluster_id)

                de_results_prot<-read.csv(de_prot_file) # read the file 
                de_results_prot_sig<-de_results_prot[de_results_prot[, metric_p]<T_p,] # take only the union of all at the end 
                mofa_features<- unique(top_proteins$feature)


               if (sig_only){

                        mofa_features <- intersect(mofa_features, de_sig_all_top)

               }
                de_results_prot_top<-de_results_prot[match(mofa_features, de_results_prot$X),]




                # TODO: print only significant 
                #print(paste('SIG in f',de_results_prot_sig$X %in% top_proteins$feature))


                #if (sig_only){
                #        # filter the ones that are de 
#
                 #       de_results_prot_top<-de_results_prot[match( na.omit(unique(de_sig_all_top)),de_results_prot$X),]

                #}
            

                # get also the pvalue 
                de_all[[cluster_id]]<-data.frame(de_results_prot_top[, metric_p])

                head(de_results_prot_top)
             


        #        print(length(de_all[[cluster_id]]))
        }

        names(de_all)
        de_all
        # TODO: add top prot
        names(de_all)<-paste0(VISIT_COMP,'_',c(1:length(clust_ids)))
        head(de_all[[4]])
        

        all_clusts_proteins_logFC<-do.call(cbind,de_all )
        head(all_clusts_proteins_logFC)

       
        all_clusts_proteins_logFC<-data.frame(all_clusts_proteins_logFC)
        colnames(all_clusts_proteins_logFC) = names(de_all)


           # all_clusts_proteins_logFC_unique<-all_clusts_proteins_logFC[!duplicated(de_results_prot_top$X),]


        all_clusts_proteins_logFC$feature<-de_results_prot_top$X



        return(all_clusts_proteins_logFC)

}
#colnames(de_results_prot)

times<-c('BL', 'V04' ,'V06', 'V08')

#colnames(results_de)

all_clusts_proteins_logFC_all_times<-lapply( times, get_de_proteins_per_tp, sig_only=sig_only,metric=metric_p_lfc,  de_sig_all_top = de_sig_all_top)

colnames(all_clusts_proteins_logFC_all_times[[1]])

all_clusts_proteins_pval_all_times<-lapply( times, get_de_proteins_per_tp,sig_only=sig_only,metric=metric_p, de_sig_all_top = de_sig_all_top )
head(all_clusts_proteins_pval_all_times)



all_clusts_times_logFC_df<-do.call(cbind, all_clusts_proteins_logFC_all_times )
all_clusts_times_pval_df<-do.call(cbind, all_clusts_proteins_pval_all_times )
dim(all_clusts_times_pval_df)
dim(all_clusts_times_logFC_df)


### REMOVE duplicated or NA gene names
valid_inds<-!is.na(all_clusts_times_logFC_df$feature ) & !duplicated(all_clusts_times_logFC_df$feature )


all_clusts_times_logFC_df<-all_clusts_times_logFC_df[valid_inds,]
all_clusts_times_pval_df<-all_clusts_times_pval_df[valid_inds,]


rownames(all_clusts_times_pval_df)<-all_clusts_times_pval_df$feature;
all_clusts_times_pval_df[colnames(all_clusts_times_pval_df) == 'feature']<-NULL

rownames(all_clusts_times_logFC_df)<-all_clusts_times_logFC_df$feature; 
all_clusts_times_logFC_df[colnames(all_clusts_times_logFC_df) == 'feature']<-NULL


x = all_clusts_times_logFC_df


all_clusts_times_logFC_df$BL_2


colnames(all_clusts_times_logFC_df)

# add factor annotation 
#top_proteins
# heatmap 
row_an<-as.factor(top_proteins$Factor[match(rownames(all_clusts_times_logFC_df),top_proteins$feature)])
row_ha<-rowAnnotation(factor=row_an)
cluster_cols=FALSE


all_clusts_times_pval_df1<-all_clusts_times_pval_df
all_clusts_times_pval_df1[all_clusts_times_pval_df<T_p]='*'
all_clusts_times_pval_df1[all_clusts_times_pval_df<(T_p/10)]='**'
all_clusts_times_pval_df1[all_clusts_times_pval_df>T_p]=''

all_clusts_times_pval_df1




gene_symbols_all<-convert_to_gene_symbol(rownames( all_clusts_times_logFC_df), view=view)

gene_symbols_all

rownames(all_clusts_times_logFC_df)<-gene_symbols_all

nf<-dim(all_clusts_times_logFC_df)[1]
nf
#cluster_params_dir
cluster_params_dir
outdir_s_p_all_vis <- paste0(cluster_params_dir, '/de_c0/')


dir.create(outdir_s_p_all_vis,'all_time/', recursive=TRUE)

all_clusts_times_logFC_df[is.na(all_clusts_times_logFC_df)]<-0
all_clusts_times_pval_df1[is.na(all_clusts_times_pval_df1)]<-''



xminxmax<-get_limits(apply(all_clusts_times_logFC_df,2, as.numeric))
all_clusts_times_logFC_df
xminxmax
col_fun = colorRamp2(c(xminxmax[1], 0, xminxmax[2]), c("blue", "white", "red"))

tissue_s<-unlist(strsplit(tissue,''))
tissue_s

prot_de_mode

column_title = paste(view_s, ',',  metric_p, '<', T_p, ',',  'DE only:', sig_only )


hname<-paste0(outdir_s_p_all_vis,'all_time/',view, tissue_s[1], '_', prot_de_mode,'_c',as.numeric(cluster_cols),'_tp_', length(times), '_',top_fr,'_s',
                 as.numeric(sig_only), 'p_', metric_p,T_p)
print(hname)
height=1+log(nf)

cluster_cols=FALSE



jpeg(paste0(hname,'_hm_log2FC.jpeg'),  res=200, width=5, height=1+log(nf), units='in')

print(paste(rownames(as.matrix(all_clusts_times_logFC_df))))
 cm<-ComplexHeatmap::pheatmap(as.matrix(all_clusts_times_logFC_df), 
 column_split = rep(clusters_names_h, length(times)), 
  col = col_fun, 
  cluster_cols = cluster_cols,
  right_annotation=row_ha,
  display_numbers = as.matrix(all_clusts_times_pval_df1), 
  

  )


 cm_all[[view]] =cm

draw(cm, column_title=column_title, 
  column_title_gp = gpar(fontsize = 10))
dev.off()
graphics.off()



#print(rownames(as.matrix(all_clusts_times_logFC_df)))
#tname<-paste0(outdir_s_p,'../all_time/',tissue,'_cc_',as.numeric(cluster_cols),'_tp_', length(times), '_',top_fr,'prot.csv')
#print(tname)
#write.csv(rownames(as.matrix(all_clusts_times_logFC_df)), tname)
print('cluster1')
View(data.frame(c(rownames(all_clusts_times_pval_df1)[which(all_clusts_times_pval_df1[, 'V08_1']!='')])))
View(c(rownames(all_clusts_times_pval_df1)[which(all_clusts_times_pval_df1[, 'V08_1']!='')]))

all_clusts_times_pval_df1
feats<-c(rownames(all_clusts_times_pval_df1)[which(all_clusts_times_pval_df1[, 'V08_1']!='')])
feats
mirs_to_enrich<-run_enrich_mirnas(feats, test_type = 'ORA',top_mirs_ora = length(feats))
colnames(mirs_to_enrich)
mirs_to_enrich$`P-adjusted`
write.csv(mirs_to_enrich_ord<-mirs_to_enrich %>% arrange(`P-adjusted`), paste0(hname, 'enrichment.csv'))
mirs_to_enrich_ord$Subcategory[1:5]
#mirna_enrich_res_postprocessing



































