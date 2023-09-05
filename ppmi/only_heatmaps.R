

cors_sig
p<-plot_data_heatmap(MOFAobject_gs, 
                     view = views[i], 
                     factor =  ii,  
                     features = nfs,
                     denoise = TRUE,
                     cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                     show_rownames = TRUE, show_colnames = TRUE,
                     scale = "row",
                     main=main_t
)
dev.off()
show(p)
     

graphics.off()
for (i in seq(1,vps)){
  for (ii in seq(1,fps)){

      cluster_rows=TRUE;cluster_cols=TRUE
      
      

      ###### Heatmaps 
      nfs=20
      #jpeg(paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs, '_cr_',cluster_rows, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
      # Plot heatmaps for each factor only for miRNA 
      
      var_captured<-round(vars_by_factor[ii,i], digits=2)
      main_t<-paste0('Factor ', ii, ', Variance = ',var_captured, '%')
                
      ns<-dim(MOFAobject@samples_metadata)[1]
      cor_T<-2; cor_p_T<-0.1
      abs(cors_pearson)>0.15
      
    dim(cors_pearson)
    dim(cors)
      rel_cors<-cors[ii,][cors[ii,]>cor_T &  cors_pearson[ii,]>cor_p_T ]
      rel_cors
      
      cors_sig=names(which(cors[ii,]>cor_T))
      FT=0
      if (length(cors_sig)==0){
        cors_sig=c()
        
      } else if (length(cors_sig)>5){
        FT=5
        # rel_cors_ordered<-rel_cors[order(-rel_cors)][1:7]
        rel_cors_ordered<-rel_cors[order(-rel_cors)]
        
        cors_sig<-names(rel_cors_ordered)
      }
      cors_sig
      exclude_vars= c('LAST_UPDATE_M4', 'INFODT_M4', 'NTEXAMTM', 'REC_ID_moca', 'REC_ID_st', 'OFFEXAMTM', 'OFFEXAMDT', 'OFFPDMEDT', 
                      'INFO')
      cors_sig<-cors_sig[!(cors_sig %in% exclude_vars)]; cors_sig
      cors_sig<-cors_sig[!grepl( 'LAST_UPDATE|INFO_DT|TM|DT|ORIG_ENTRY|DATE|PAG', cors_sig)]

      plot_heatmap_flag=TRUE
      #MOFAobject_gs@samples_metadata[cors_sig][is.na(MOFAobject_gs@samples_metadata[cors_sig])]<-10^-6v
      MOFAobject_gs@samples_metadata[,cors_sig]
      MOFAobject_gs@samples_metadata[,cors_sig]
      
      cors_sig_non_na<-names(which(!apply(is.na(MOFAobject_gs@samples_metadata[,cors_sig]),2,any )))
      if(length(cors_sig_non_na)==0){
        cors_sig_non_na=c()
      }
        hname<-paste0(outdir, 'heatmap/heatmap_',ii,'_',views[i],'_', 'nfs_', nfs,'_cr_', cluster_rows, res, '_cor_', cor_T, 'FT_', FT, '.jpeg')
        
        MOFAobject_gs@samples_metadata[cors_sig_non_na]
        p<-plot_data_heatmap(MOFAobject_gs, 
                             view = views[i], 
                             factor =  ii,  
                             features = nfs,
                             denoise = TRUE,
                             cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                             show_rownames = TRUE, show_colnames = TRUE,
                             scale = "row",
                             annotation_samples = cors_sig_non_na,
                             main=main_t
                             
                             
        )
        #ggsave(hname, plot=p,height=nfs/2, width=(ns+as.numeric(length(cors_sig_non_na) )) )
        ggsave(hname, plot=p,height=nfs/2, width=ns/50, dpi=250) 
        
      
  }
    # top weights
    # concat all 
    
    
    
  }
  
  

