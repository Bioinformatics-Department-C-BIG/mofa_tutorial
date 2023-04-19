
#### RNA comparison 

rna_p<-paste0('ppmi/plots/single/rnas_', VISIT, '_0.1_100_coh_1-2_AGE_AT_VISIT+SEX+COHORT/gseGO_BP_1_0log2pval.csv')
rna_paths<-read.csv(rna_p)
rna_paths_cut<-rna_paths[rna_paths$p.adjust<Padj_T_paths,]
dim(rna_paths_cut)




mofa_p<-paste0('ppmi/plots/p_',VISIT,'_Plasma_0.9_T_1-2vsn_TNA_0.8g_0.5_100_m_0.75_10_8_coh_1-2_V08_TRUE/enrichment/GO_BP_0.05_enrichment_positive_pvals_no_f.csv')
mofa_paths<-read.csv(mofa_p)

mofa_paths_list<-mofa_paths[,3]
mofa_paths_list<-tolower(gsub('_', ' ',gsub('GOBP_', '', mofa_paths_list)))
mofa_paths_list
# todo: convert everything to caps to comapre? ???


#### COMPARE ALSO WITH MULTI!!! 

if (grepl('GO', mir_paths[1])){
  if (startsWith(mir_paths[1], 'GO')){
    mir_paths<-sub(".*? ", "", mir_paths)
    
  }else{
    mir_paths<-gsub("\\ GO.*", "", mir_paths)
    
  }
  
}

mir_paths[!(mir_paths %in% common_pathways)]
common_pathways<-intersect(rna_paths_cut$Description ,mir_paths)

length(common_pathways)


listInput=list(RNA=rna_paths_cut$Description, 
               miR=mir_paths, 
               mofa=mofa_paths_list)


library(RColorBrewer)
ni<-length(listInput)
myCol <- brewer.pal(ni, "Pastel2")[1:ni]

venn.diagram(listInput,   
             filename = paste0(results_file, '_bar_14_venn_diagramm.png'), 
             
             
             fill=myCol,
             output=TRUE)


dev.off()



#dev.off()
#, color=P-adjusted)

