




### Heatmaps of the log2FC values 
# Input: log2FC from deseq with controls for each gene for each time points 
# Read in deseq results from each timepoint 
# And each cluster 
# need to run after DE tutorial ppmi cluster_compare.R

#write.csv(results_de, paste0(outdir_s_p, 'results.csv'))
#VISIT_COMP='V06'
prefix='prot_'
DIFF_VAR
# TODO: run all timepoints!!!
tissue_un_mofa<-ifelse(tissue_un=='Plasma', 'plasma', 'csf')
tissue_un

#' @param tissue is a global name that adjusts for targeted or untargeted 
if (prot_de_mode=='t'){
        view=paste0('proteomics_t_', tolower(TISSUE))
        top_fr=0.05
        tissue=TISSUE
}else{
        view=paste0('proteomics_', tolower(tissue_un_mofa))
        print(view)
        tissue = tissue_un
        top_fr=0.025

}
top_fr=0.05
#fact2 = c(13,22)
top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact, view=view, top_fr=top_fr   )
#top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact2, view=view, top_fr=top_fr   )

top_proteins
top_proteins$feature<-gsub(paste0('_',view),'', top_proteins$feature)
top_proteins$feature


clust_ids<-c('1', '2', '3', '4')




#cluster_params_dir
# outdir_s_p
fact = get_factors_for_metric(DIFF_VAR)
fact
cluster_params_dir<-get_cluster_params_dir(DIFF_VAR)
cluster_params_dir




metric_p<-'adj.P.Val'
metric_p='P.Value'; T_p=0.01

de_sig_all<-list()
  for (cluster_id in clust_ids){

                outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
                de_prot_file<-paste0(outdir_s_p, prefix, tissue,'_', prot_de_mode,'_de_cl',cluster_id,  '_results.csv')
                de_results_prot<-read.csv(de_prot_file)
                
                view=paste0('proteomics_', tolower(TISSUE))
                # TODO: take the top MOFA proteins from moca 
              
                #}
                de_results_prot_sig<-de_results_prot[de_results_prot$P.Value<0.01,] # take only the union of all at the end 
                de_results_prot_sig<-de_results_prot[de_results_prot[,metric_p]<T_p,] # take only the union of all at the end 


                de_sig_all[[cluster_id]]<-de_results_prot_sig$X

  }
de_sig_all_top<-unique(unlist(de_sig_all))
de_sig_all_top
# TODO: separate to get top 

get_de_proteins_per_tp<-function(VISIT_COMP, metric_p='logFC', sig_only =FALSE){
        de_all<-list()

        print(paste(VISIT_COMP, metric_p))
        for (cluster_id in clust_ids){

                outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
                de_prot_file<-paste0(outdir_s_p, prefix, tissue,'_', prot_de_mode,'_de_cl',cluster_id,  '_results.csv')
                de_results_prot<-read.csv(de_prot_file)
                
                view=paste0('proteomics_', tolower(TISSUE))
                # TODO: take the top MOFA proteins from moca 
              
                #}
                #de_results_prot_sig<-de_results_prot[de_results_prot$P.Value<0.01,] # take only the union of all at the end 
                de_results_prot_sig<-de_results_prot[de_results_prot[, metric_p]<T_p,] # take only the union of all at the end 

                print(paste('SIG',de_results_prot_sig$X))
                de_results_prot_top<-de_results_prot[de_results_prot$X %in% top_proteins$feature,]
                # TODO: print only significant 
                print(paste('SIG in f',de_results_prot_sig$X %in% top_proteins$feature))

                if (sig_only){
#
                        de_results_prot_top<-de_results_prot[de_results_prot$X %in% de_sig_all_top,]
                }

               # de_results_prot_top<-de_results_prot[de_results_prot$X %in% de_results_prot_sig$X,]

                # get also the pvalue 
                de_all[[cluster_id]]<-de_results_prot_top[, metric_p]

        #        print(length(de_all[[cluster_id]]))
        }
        

        # TODO: add top prot
        names(de_all)<-paste0(VISIT_COMP,'_',c(1:length(clust_ids)))
        all_clusts_proteins_logFC<-do.call(cbind,de_all )
        all_clusts_proteins_logFC<-as.data.frame(all_clusts_proteins_logFC)
        rownames(all_clusts_proteins_logFC)<-de_results_prot_top$X
        colnames(all_clusts_proteins_logFC)<-names(de_all)


        return(all_clusts_proteins_logFC)

}
#colnames(de_results_prot)

times<-c('BL', 'V04' ,'V06', 'V08')


#colnames(results_de)
sig_only=TRUE
all_clusts_proteins_logFC_all_times<-lapply( times, get_de_proteins_per_tp, sig_only=sig_only)
all_clusts_proteins_logFC_all_times

all_clusts_proteins_pval_all_times<-lapply( times, get_de_proteins_per_tp,sig_only=sig_only,metric=metric_p )
all_clusts_proteins_pval_all_times



all_clusts_times_logFC_df<-do.call(cbind, all_clusts_proteins_logFC_all_times )
all_clusts_times_pval_df<-do.call(cbind, all_clusts_proteins_pval_all_times )
dim(all_clusts_times_pval_df)
dim(all_clusts_times_logFC_df)

x = all_clusts_times_logFC_df





colnames(all_clusts_times_logFC_df)
# add factor annotation 
#top_proteins
row_an<-as.factor(top_proteins$Factor[match(rownames(all_clusts_times_logFC_df),top_proteins$feature)])
row_ha<-rowAnnotation(factor=row_an)
cluster_cols=FALSE


all_clusts_times_pval_df1<-all_clusts_times_pval_df
all_clusts_times_pval_df1[all_clusts_times_pval_df<T_p]='*'
all_clusts_times_pval_df1[all_clusts_times_pval_df>T_p]=''

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
prot_de_mode
all_clusts_times_logFC_df[is.na(all_clusts_times_logFC_df)]<-0
all_clusts_times_pval_df1[is.na(all_clusts_times_pval_df1)]<-''

all_clusts_times_logFC_df

xminxmax<-get_limits(all_clusts_times_logFC_df)
xminxmax
col_fun = colorRamp2(c(xminxmax[1], 0, xminxmax[2]), c("blue", "white", "red"))

tissue_s<-unlist(strsplit(tissue,''))



hname<-paste0(outdir_s_p_all_vis,'all_time/',tissue_s[1], '_', prot_de_mode,'_cc_',as.numeric(cluster_cols),'_tp_', length(times), '_',top_fr,'_s',
                 as.numeric(sig_only), 'p_', metric_p,'_hm_log2FC.jpeg')
print(hname)
height=1+log(nf)

jpeg(hname,  res=200, width=5, height=1+log(nf), units='in')

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






































