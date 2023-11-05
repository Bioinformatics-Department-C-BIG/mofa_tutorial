


if (view=='RNA'){
  ens<-as.character(de_group_vs_control1$symbol)
  symb<-get_symbols_vector(ens)
  de_group_vs_control1$GENE_SYMBOL<-symb
  #feat_names_ens_ids<-unique(symb)

}


#### OVERLAP OF A LIST OF MIRS AND THEIR TARGETS
## TIME related rnas or DE for group 1?? 
list_of_rnas<-most_sig_over_time$GENE_SYMBOL
###
list_of_rnas<-de_group_vs_control1$GENE_SYMBOL
list_of_rnas

#gsea_results_fname<-paste0(mir_results_file,'_mieaa_res.csv' )
pvalueCutoff=1
mirs=gsub( '\\.','-', de_group_vs_control_and_time2)
mirs=gsub( '\\.','-', selected_mirs)
## NOTE : the more mirs you add the less gene targets you get back as significant....
test_type='ORA'
if (file.exists(gsea_results_fname)){
  ### Load enrichment results if available
  mieaa_all_gsea<-read.csv(gsea_results_fname, header=TRUE)
  ### TODO: Rerun with updated pvalue cutoff 
}else{
  ## otherwise run GSEA analysis 
  
  
  mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                                     mirna_type = "mature",
                                     test_type = test_type,
                                     species = 'Homo sapiens',
                                     sig_level=pvalueCutoff
  )
  
  
  write.csv(mieaa_all_gsea, gsea_results_fname, row.names = FALSE)
  
}

colnames(mieaa_all_gsea)<-make.names(  colnames(mieaa_all_gsea))
mieaa_all_gsea_sig<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('GO Biological process (miRPathDB)')) %>%
  dplyr::filter(P.adjusted<0.05)

mieaa_all_gsea_sig$Subcategory
#View(mieaa_all_gsea_sig)

mieaa_targets<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('Target genes (miRTarBase)')) %>%
  dplyr::filter(P.value<0.05)

head(mieaa_targets)
#View(mieaa_targets)


#####
intersect(list_of_rnas, mieaa_targets$Subcategory)

mieaa_targets$Subcategory






