

#library('MOFAdata')
library('MOFA2')



## TODO: if enrichment is already run then just load results
## load res.positive to be used in the next script
#res.positive$feature.sets
#res=res.positive
#es.positive$pval.adj
#res=res.negative

write_enrich<-function(res, sign_mode){
  #' 
  #'' @res res.negative result from mofa enrichment 
  #'
  #'
  #'res$pval.adj
  #' res
      colnames(res)
      head(res)
      res$pval.adj=as.data.frame(res$pval.adj)   ;
      results_enrich<-as.data.frame(sapply(res$pval.adj, as.numeric));
      all_fs<-colnames(results_enrich)
      results_enrich$Description=rownames(res$pval.adj); 

      # sapply(all_fs,pcgse_dot_by_factor, results_enrich=results_enrich ) # do not run dotplot
      



      all_fs_merged2<-reshape::melt(results_enrich, id.vars = 'Description');
      
      res$pval=as.data.frame(res$pval) 
      res_enrich_pval=as.data.frame(sapply(res$pval, as.numeric))
      res_enrich_pval$Description=rownames(res$pval);
      
      all_fs_merged2_pval<-reshape::melt(res_enrich_pval, id.vars = 'Description');
      all_fs_merged2_pval2=cbind(all_fs_merged2, all_fs_merged2_pval$value)
      #rownames(all_fs_merged2_pval2)=rownames(all_fs_merged2_pval);
      
      #all_fs_merged2_pval2<-merge(all_fs_merged2,all_fs_merged2_pval , by=c('variable', 'value'))

      rownames(all_fs_merged2_pval2)
      all_fs_merged2_pval2<-all_fs_merged2_pval2[with(all_fs_merged2_pval2, order(variable, value)),]# order 
     colnames(all_fs_merged2_pval2)<-c('Description', 'Factor', 'p.adjust', 'pvalue')
      #colnames(all_fs_merged2_pval2)<-c('Factor', 'p.adjust', 'pvalue')
      
      neg_file<-paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                       mode, '_enrichment', sign_mode)
      
      write.csv(format(all_fs_merged2_pval2, digits=3),paste0(neg_file,  '.csv' ), row.names = TRUE)
      all_fs_merged2_pval2_sig=all_fs_merged2_pval2[ all_fs_merged2_pval2$p.adjust<T,]
      write.csv(format(all_fs_merged2_pval2_sig, digits=3),paste0(neg_file, '_', T,  '.csv' ))
      return(all_fs_merged2_pval2_sig)
      
}

# plot: dot plot
pcgse_dot_by_factor<-function(factor, results_enrich){

      barpl_input<-results_enrich[ c('Description',factor )]
      colnames(barpl_input)<-c('Description','p.adjust' )
     barpl_input<-barpl_input[barpl_input$p.adjust<0.05,]

      barpl_input_top<-barpl_input[order(barpl_input[,2], decreasing=FALSE)[1:20],]

    barpl_input_top$x = factor
    barpl_input_top$log10=-log10(barpl_input_top$p.adjust)
    barpl_input_top$log10
    ggplot(data = barpl_input_top, aes( x=x,y = Description, 
                            color = `p.adjust`, size=-log10(p.adjust))) + 
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      theme_bw() + 
      ylab("") + 
      xlab("") + 
      ggtitle("GO enrichment analysis")


    dir.create(paste0(outdir, '/enrichment/pcgse/',mode, '/'))
    ggsave(paste0(outdir, '/enrichment/pcgse/',mode, '/', subcategory_s,'_', sign_mode, '_dp_', factor, '.png'), width=7, height=5)
    }


subcategory<- 'CP:KEGG'
subcategory<- 'CP:KEGG'
subcategory<- 'GO:MF'
subcategory<- 'GO:BP'
dir.create(paste0(outdir, '/enrichment/'))
#for (subcategory in c('GO:BP' ,'CP:KEGG')){

mode='proteomics'
mode='proteomics_csf'


assay(mofa_multi[,,3])



mode='proteomics_t_csf'

   
features_names(MOFAobject)$RNA
features_names(MOFAobject)$RNA<-sapply(features_names(MOFAobject)$RNA, 
                                         function(x) {stringr::str_remove(x, '\\..*')})

MOFAobject_enr<-MOFAobject
features_names(MOFAobject_enr)$proteomics_csf<-gsub('_proteomics_csf','',features_names(MOFAobject_enr)$proteomics_csf)
features_names(MOFAobject_enr)$proteomics_csf
features_names(MOFAobject_enr)$proteomics_plasma<-gsub('_proteomics_plasma','',features_names(MOFAobject_enr)$proteomics_plasma)
features_names(MOFAobject_enr)$proteomics_t_plasma<-gsub('_proteomics_t_plasma','',features_names(MOFAobject_enr)$proteomics_t_plasma)
features_names(MOFAobject_enr)$proteomics_t_csf<-gsub('_proteomics_t_csf','',features_names(MOFAobject_enr)$proteomics_t_csf)
#features_names(MOFAobject_enr)$proteomics_t_csf
library('org.Hs.eg.db')
library('AnnotationDbi')

mode='proteomics_t_csf'
mode = 'proteomics_csf'
mode = 'proteomics_t_plasma'


mode = 'proteomics_csf'

mode='RNA'

my_protein_ids = features_names(MOFAobject_enr)[[mode]]

head(my_protein_ids)
#my_protein_mapping<-AnnotationDbi::select(org.Hs.eg.db, my_protein_ids,  "SYMBOL","UNIPROT")
#my_protein_mapping<-my_protein_mapping[!duplicated(my_protein_mapping$UNIPROT),]
#my_protein_mapping<-my_protein_mapping[match(features_names(MOFAobject_enr)[[mode]], my_protein_mapping$UNIPROT),]
#features_names(MOFAobject_enr)[[mode]]<-my_protein_mapping$SYMBOL


force_pcgse = TRUE
get_feature_set_uniprot<-function(gs_original){
  #' map colnames of feature.sets with UNIRPOT SYMBOL
  #' @param gs_original: the feature set with gene symbols 
  #' 
      gs_uniprot<-gs_original
      head(colnames(gs_original))
    # convert the feature.sets matrix rownames to uniprot
    #colnames(gs_uniprot)<-
    ug_mapping<-AnnotationDbi::select(org.Hs.eg.db, colnames(gs_original),"UNIPROT", "SYMBOL")
    # filter by what is found
    na.omit(ug_mapping$UNIPROT)
    ug_mapping_found<-ug_mapping[!is.na(ug_mapping$UNIPROT),]
    ug_mapping_found$UNIPROT
    colnames_found<-colnames(gs_uniprot)[colnames(gs_uniprot) %in%  ug_mapping_found$SYMBOL]
    gs_uniprot_filt<-gs_uniprot[,colnames_found]
    length(ug_mapping_found$UNIPROT)
    colnames(gs_uniprot_filt)<-ug_mapping_found$UNIPROT[match(colnames(gs_uniprot_filt),ug_mapping_found$SYMBOL )]
    colnames(gs_uniprot_filt)
    return(gs_uniprot_filt)
}







#rownames(results_de)
grepl('proteomics', mode)
# 'GO:MF'
for (subcategory in c('GO:BP' )){

  if (grepl('proteomics', mode)){
      gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')
      gs_original<-as.matrix(read.csv(gs_file, header=1, row.names=1))
      gs<-gs_original
    if (mode=='proteomics_csf'|| mode == 'proteomics_plasma'){
      # untargeted are in uniprot..
      gs<-get_feature_set_uniprot(gs_original)
    }

  }else{
    #output_files
      gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
      #gs_file
      gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))
   
  }




  sign_mode='negative'
  subcategory_s<-gsub('\\:', '_', subcategory)
  enrich_res_file_neg<-paste0(outdir,'/enrichment/pcgse/' ,subcategory_s, '_', T, mode, '_enrichment_', 'negative' )
  enrich_res_file_pos<-paste0(outdir,'/enrichment/pcgse/' ,subcategory_s, '_', T, mode, '_enrichment_', 'positive' )
  

  if (!force_pcgse & file.exists(enrich_res_file_neg)){
        res.negative=loadRDS(enrich_res_file_neg)
        res.positive=loadRDS(enrich_res_file_pos)
    
  }else{
    #res.negative
          # GSEA on negative weights, with default options
          res.negative <- run_enrichment(MOFAobject_enr, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "negative"
          )
          
          res.positive <- run_enrichment(MOFAobject_enr, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "positive"
          )
          
          sign_mode='negative'
          res_negative_df<-write_enrich(res.negative, sign_mode=sign_mode)
          res_negative_df
          saveRDS(res.negative, paste0(outdir,'/enrichment/pcgse/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'negative' ))
          
          sign_mode='positive'

          res_positive_df<-write_enrich(res.positive, sign_mode=sign_mode)
          saveRDS(res.positive, paste0(outdir,'/enrichment/pcgse/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'positive' ))
          res_negative_df
          res_merged<-merge(res_negative_df, res_positive_df, suffixes=c('_n', '_p'),
            by=c('Description','Factor' ))
          #res_merged$pvalue_min<-c(rowMins(as.matrix(res_merged[,c('pvalue_n','pvalue_p')])))
          
          write.csv(res_merged, paste0(outdir,'/enrichment/pcgse/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment.csv' ), row.names=FALSE)
  }
  
  
  
  

    ## TODO: create a function to do for both positive and negative 
    #
    T=0.05
  

  ##### which factor is related to parkinsons disease in KEGG
  ### PROBLEM: this is based on RNA only!!! 
  alpha=0.05
dir.create(paste0(outdir, '/enrichment/pcgse/', mode,'/'), recursive=TRUE)

## After having run with uniprot, convert the results to symbol before plotting ####




  convert_to_gene_symbol<-function( x, view)(
    #'view = proteomics csf/ proteomics targeted/rna 
    #' 
    #' 
  

    if (view == 'RNA'){
      x = get_symbols_vector(x)
    } else if (view %in% c('proteomics_csf', 'proteomics_plasma')){
      x = get_symbol_from_uniprot(x)$SYMBOL

    }
  
  
  )

  rownames(res.negative$feature.statistics) = convert_to_gene_symbol(rownames(res.negative$feature.statistics), view=mode)
  colnames(res.negative$feature.sets)<-convert_to_gene_symbol(colnames(res.negative$feature.sets),  view=mode)


  rownames(res.positive$feature.statistics)=convert_to_gene_symbol(rownames(res.positive$feature.statistics),  view=mode)
  colnames(res.positive$feature.sets)<-convert_to_gene_symbol(colnames(res.positive$feature.sets),  view=mode)





sapply(1:N_FACTORS, function(factor){
  tryCatch({

  plot_enrichment_detailed(res.negative, factor, 
  alpha = alpha, 
  max.genes=6
  )
  #graphics.off()
  ggsave(paste0(outdir, '/enrichment/pcgse/', mode,'/',subcategory_s, '_detailed_neg', '_', factor, '.png'),
  width=6, height=4)

  },
  error = function(e) {an.error.occured <<- TRUE}
  )
tryCatch({
 

  plot_enrichment_detailed(res.positive, factor, 
    max.genes=6,

    alpha = alpha)
  ggsave(paste0(outdir, '/enrichment/pcgse/', mode, '/', subcategory_s, '_detailed_pos', '_', factor, '.png'),
  width=6, height=4)

  },
  error = function(e) {an.error.occured <<- TRUE}
  )


}



)

  
}




sign_mode='negative'
subcategory<- 'GO:BP'
T=0.05


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


#res.negative %>% 
#  dplyr::filter(pval.adj<0.05)
plot_enrichment_heatmap(res.negative, 
                        alpha=0.00000000004, cap=0.00000000005 
)

dev.off()




###### turn to enrichment result to plot 

# MERGE RESULTS CSV

all_modes<-c('proteomics_csf', 'proteomics_plasma', 'RNA')
 all_enrichments =  list()

i=1
for (mode in all_modes){
  for (sign_mode in c('positive', 'negative')){
      
       neg_file<-paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                       mode, '_enrichment', sign_mode,'_', T,  '.csv')

        results_enrich_tmp<-read.csv(neg_file)
        results_enrich_tmp$sign_mode = sign_mode
        results_enrich_tmp$mode = mode
        all_enrichments[[i]] = results_enrich_tmp
        i=i+1
#        print(head(results_enrich_tmp))


  }


}

 
all_enrichments_bind<-do.call(bind_rows,all_enrichments)

write.csv(all_enrichments_bind, paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                        '_enrichment','_', T,'all_modes',   '.csv'))






















