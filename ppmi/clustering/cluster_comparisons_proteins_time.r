




### Heatmaps of the log2FC values 
# Input: log2FC from deseq with controls for each gene for each time points 
# Read in deseq results from each timepoint 
# And each cluster 

#write.csv(results_de, paste0(outdir_s_p, 'results.csv'))
VISIT_COMP='V06'
TISSUE
# TODO: run all timepoints!!!
view=paste0('proteomics_', tolower(TISSUE))
top_proteins<-concatenate_top_features(MOFAobject, factors_all=fact, view=view, top_fr=0.01   )
top_proteins$feature<-gsub(paste0('_',view),'', top_proteins$feature)
top_proteins$feature


get_de_proteins_per_tp<-function(VISIT_COMP){
        de_all<-list()


        for (cluster_id in c(1:3)){

            outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
            de_prot_file<-paste0(outdir_s_p, prefix,TISSUE,'_de_cl',cluster_id,  '_results.csv')
            de_results_prot<-read.csv(de_prot_file)
            
            cluster_id
                view=paste0('proteomics_', tolower(TISSUE))
                # TODO: take the top MOFA proteins from moca 
                paste0('_',view)

                top_proteins
                de_results_prot$X
                top_proteins$feature
                de_results_prot_top<-de_results_prot[de_results_prot$X %in% top_proteins$feature,]
                

                de_all[[cluster_id]]<-de_results_prot_top$logFC
                print(length(de_all[[cluster_id]]))
        }
        


        names(de_all)<-paste0(VISIT_COMP,'_',c(1:3))
        all_clusts_proteins_logFC<-do.call(cbind,de_all )
        all_clusts_proteins_logFC
        all_clusts_proteins_logFC<-as.data.frame(all_clusts_proteins_logFC)
        rownames(all_clusts_proteins_logFC)<-de_results_prot_top$X
        de_results_prot_top$X
        colnames(all_clusts_proteins_logFC)<-names(de_all)
        print(all_clusts_proteins_logFC)
        
        return(all_clusts_proteins_logFC)

}



all_clusts_proteins_logFC_all_times<-lapply(c('BL', 'V04', 'V06', 'V08'), get_de_proteins_per_tp)
all_clusts_proteins_logFC_all_times_df<-do.call(cbind, all_clusts_proteins_logFC_all_times )

VISIT_COMP
rep(1:3, 4)

colnames(all_clusts_proteins_logFC_all_times_df)
jpeg(paste0(outdir_s_p,'/',TISSUE, '_heatmap_log2FC.jpeg'),  res=200, width=10, height=10, units='in')
ComplexHeatmap::pheatmap(as.matrix(all_clusts_proteins_logFC_all_times_df), 
 column_split = rep(1:3, 4), )
dev.off()













