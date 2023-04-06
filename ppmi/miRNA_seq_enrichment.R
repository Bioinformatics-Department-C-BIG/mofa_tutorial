
#install.packages('rbioapi')
library(rbioapi)
order_by_metric<-'abslog2pval'
if (VISIT=='V08'){
  padj_T=0.01
  log2fol_T=0.1
}else{
  padj_T=0.05
  log2fol_T=0.1
}
  


library(enrichplot)


gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T, log2fol_T )
top_n=length(gene_list);
#top_n=200
top_n

gene_list_cut<-gene_list[1:top_n]
length(gene_list); length(gene_list_cut)
mirs=names(gene_list_cut)

mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                              mirna_type = "mature",
                              test_type = "GSEA",
                              species = 'Homo sapiens'
                              )

table(mieaa_all_gsea$Category)
Category<-'Annotation (Gene Ontology)';
#Category<-'GO Biological process (miRPathDB)';
#Category<-'Gene Ontology (miRWalk)';


mieaa_gsea_1<-mieaa_all_gsea[mieaa_all_gsea$Category==Category,]

Padj_T_paths=0.01
mieaa_gsea_1_cut<-mieaa_gsea_1[mieaa_gsea_1$`P-adjusted`<Padj_T_paths, ]
mir_paths<-mieaa_gsea_1_cut[,c(2)]

results_file<-paste0(outdir_s, '/mirs_enrich_', '_', padj_T, '_',  log2fol_T, '_',  order_by_metric, '_',top_n, '_', Category)

write.csv(mieaa_all_gsea, paste0(results_file, '_', '.csv' ))



##### plot results 
##
#DOSE::enrichResult
#  EnrichResult(mieaa_all_gsea)

  

df=mieaa_gsea_1
df$`P-adjusted`<-as.numeric(df$`P-adjusted`)
df$padj<-as.numeric(df$`P-adjusted`)

df$Observed<-as.numeric(df$Observed)

df_ord<-df[order(df$padj),]

df_ord<-df_ord[1:70,]
df_ord$padj

ggplot(df_ord, aes(x=reorder(Subcategory, padj), y=Observed, fill=padj))+
  geom_bar(position='dodge', stat='identity')+
  coord_flip()
ggsave(paste0(results_file, '_bar',  '.png'), height = 7, width=8)




############

#### RNA comparison 

rna_p<-paste0('ppmi/plots/single/rnas_', VISIT, '_0.1_100_coh_1-2_AGE_AT_VISIT+SEX+COHORT/gseGO_BP_1_0log2pval.csv')
rna_paths<-read.csv(rna_p)


rna_paths_cut<-rna_paths[rna_paths$p.adjust<Padj_T_paths,]
dim(rna_paths_cut)


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
               miR=mir_paths)


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

