




### Heatmaps of the log2FC values 
# Input: log2FC from deseq with controls for each gene for each time points 
# Read in deseq results from each timepoint 
# And each cluster 
# need to run after DE tutorial ppmi cluster_compare.R

#write.csv(results_de, paste0(outdir_s_p, 'results.csv'))
#VISIT_COMP='V06'
prefix='prot_'
# TODO: run all timepoints!!!
tissue_un_mofa<-ifelse(tissue_un=='Plasma', 'plasma', 'csf')
tissue_un

#' @param tissue is a global name that adjusts for targeted or untargeted 
if (prot_de_mode=='t'){
        view=paste0('proteomics_t_', tolower(TISSUE))
        top_fr=0.01
        tissue=TISSUE
}else{
        view=paste0('proteomics_', tolower(tissue_un_mofa))
        tissue = tissue_un
        top_fr=0.02

}

top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact, view=view, top_fr=top_fr   )
top_proteins
top_proteins$feature<-gsub(paste0('_',view),'', top_proteins$feature)
top_proteins$feature


clust_ids<-c('1', '2', '3')




#cluster_params_dir
# outdir_s_p
fact = get_factors_for_metric(DIFF_VAR)
fact
cluster_params_dir<-get_cluster_params_dir(DIFF_VAR)
cluster_params_dir


# TODO: top proteins in 
get_de_proteins_per_tp<-function(VISIT_COMP, metric='logFC'){
        de_all<-list()


        for (cluster_id in clust_ids){

                outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
                de_prot_file<-paste0(outdir_s_p, prefix,tissue,'_de_cl',cluster_id,  '_results.csv')
                de_results_prot<-read.csv(de_prot_file)
                
                view=paste0('proteomics_', tolower(TISSUE))
                # TODO: take the top MOFA proteins from moca 
               
                de_results_prot_sig<-de_results_prot[de_results_prot$adj.P.Val<0.05,] # take only the union of all at the end 
                
                de_results_prot_top<-de_results_prot[de_results_prot$X %in% top_proteins$feature,]
                
                # get also the pvalue 
                de_all[[cluster_id]]<-de_results_prot_top[, metric]

                print(length(de_all[[cluster_id]]))
        }
        


        names(de_all)<-paste0(VISIT_COMP,'_',c(1:3))
        all_clusts_proteins_logFC<-do.call(cbind,de_all )
        all_clusts_proteins_logFC
        all_clusts_proteins_logFC<-as.data.frame(all_clusts_proteins_logFC)
        rownames(all_clusts_proteins_logFC)<-de_results_prot_top$X
        de_results_prot_top$X
        colnames(all_clusts_proteins_logFC)<-names(de_all)
        
        return(all_clusts_proteins_logFC)

}
#colnames(de_results_prot)

times<-c('BL', 'V04' ,'V06', 'V08')


metric_p<-'adj.P.Val'
metric_p=='P.Value'
#colnames(results_de)
all_clusts_proteins_logFC_all_times<-lapply( times, get_de_proteins_per_tp)
all_clusts_proteins_pval_all_times<-lapply( times, get_de_proteins_per_tp,metric=metric_p )


all_clusts_times_logFC_df<-do.call(cbind, all_clusts_proteins_logFC_all_times )
all_clusts_times_pval_df<-do.call(cbind, all_clusts_proteins_pval_all_times )
dim(all_clusts_times_pval_df)
dim(all_clusts_times_logFC_df)

x = all_clusts_times_logFC_df
get_limits<-function(x){
        xmax<-max(x)
        xmin<-min(x)
        max_all<-max(c(abs(xmax), abs(xmin)))
        return(c(-max_all, max_all))
}

xminxmax<-get_limits(all_clusts_times_logFC_df)
col_fun = colorRamp2(c(xminxmax[1], 0, xminxmax[2]), c("blue", "white", "red"))


colnames(all_clusts_times_logFC_df)
# add factor annotation 

row_an<-as.factor(top_proteins$Factor[match(rownames(all_clusts_times_logFC_df),top_proteins$feature)])
row_ha<-rowAnnotation(factor=row_an)
cluster_cols=FALSE


all_clusts_times_pval_df1<-all_clusts_times_pval_df
all_clusts_times_pval_df1[all_clusts_times_pval_df<0.05]='*'
all_clusts_times_pval_df1[all_clusts_times_pval_df>0.05]=''

all_clusts_times_pval_df1



 #uniprot_ids<-rownames( all_clusts_times_logFC_df)
 #gene_symbols<-get_symbol_from_uniprot(uniprot_ids)


 
   if (prot_de_mode == 'u'){
    uniprot_ids<-rownames( all_clusts_times_logFC_df)
    gene_symbols_all<-get_symbol_from_uniprot(uniprot_ids)
    gene_symbols<-gene_symbols_all$SYMBOL
   }else{
    gene_symbols<-rownames( all_clusts_times_logFC_df)
   }

rownames( all_clusts_times_logFC_df)<-gene_symbols

nf<-dim(all_clusts_times_logFC_df)[1]
nf
#cluster_params_dir
cluster_params_dir
outdir_s_p_all_vis <- paste0(cluster_params_dir, '/de_c0/')


dir.create(outdir_s_p_all_vis,'all_time/', recursive=TRUE)
hname<-paste0(outdir_s_p_all_vis,'all_time/',tissue, '_', prot_de_mode,'_cc_',as.numeric(cluster_cols),'_tp_', length(times), '_',top_fr, '_heatmap_log2FC.jpeg')

jpeg(hname,  res=200, width=6, height=3+log(nf), units='in')

cm<-ComplexHeatmap::pheatmap(as.matrix(all_clusts_times_logFC_df), 
 column_split = rep(1:length(clust_ids), length(times)), 
  col = col_fun, 
  cluster_cols = cluster_cols,
  right_annotation=row_ha,
  display_numbers = as.matrix(all_clusts_times_pval_df1)

  )
 show(cm)
dev.off()
graphics.off()



#print(rownames(as.matrix(all_clusts_times_logFC_df)))
#tname<-paste0(outdir_s_p,'../all_time/',tissue,'_cc_',as.numeric(cluster_cols),'_tp_', length(times), '_',top_fr,'prot.csv')
#print(tname)
#write.csv(rownames(as.matrix(all_clusts_times_logFC_df)), tname)





































