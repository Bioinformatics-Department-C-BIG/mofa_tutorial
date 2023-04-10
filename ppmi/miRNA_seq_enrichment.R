
#install.packages('VennDiagram')
library(rbioapi)
library('VennDiagram')
library(enrichplot)



order_by_metric<-'abslog2pval'
order_by_metric<-'abslog2pval'


if (VISIT=='V08'){
  padj_T=0.01
  log2fol_T=0.1
  padj_T=0.05
  log2fol_T=0
}else{
  padj_T=0.05
  log2fol_T=0.1
}
  



### run the enrichment if it has not been ran yet!! 

Padj_T_paths=0.01

deseq2ResDF = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)

gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T, log2fol_T )
top_n=length(gene_list);
#top_n=200

gene_list_cut<-gene_list[1:top_n]
length(gene_list); length(gene_list_cut)
mirs=names(gene_list_cut)
mirs
mir_results_file<-paste0(outdir_enrich, '/mirs_enrich_', '_', padj_T, '_',  log2fol_T, '_',  order_by_metric, '_',top_n)


gsea_results_fname<-paste0(mir_results_file,'_mieaa_res.csv' )
  
  

if (file.exists(gsea_results_fname)){
  ### Load enrichment results if available
  mieaa_all_gsea<-read.csv(gsea_results_fname, header=TRUE)
  
}else{
  ## otherwise run GSEA analysis 
  mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                                     mirna_type = "mature",
                                     test_type = "GSEA",
                                     species = 'Homo sapiens'
  )
  ### write to file to load next time 
  write.csv(mieaa_all_gsea, gsea_results_fname, row.names = FALSE)
  
}


## remove . from names 
colnames(mieaa_all_gsea)<-gsub('-','.', colnames(mieaa_all_gsea))




#table(mieaa_all_gsea$Category)
Category<-'Annotation (Gene Ontology)';
Category<-'GO Biological process (miRPathDB)';
#Category<-'Gene Ontology (miRWalk)';


mieaa_gsea_1<-mieaa_all_gsea[mieaa_all_gsea$Category==Category,]

mieaa_gsea_1_cut<-mieaa_gsea_1[mieaa_gsea_1$P.adjusted<Padj_T_paths, ]
mir_paths<-mieaa_gsea_1_cut[,c(2)]


results_df<-paste0(mir_results_file, '_', Category)

write.csv(mieaa_all_gsea, paste0(mir_results_file, '_', '.csv' ))



### ANOTHER WAY to obtain targets 
##
##

#BiocManager::install("targetscan.Hs.eg.db")
### first retrieve all targets that mieaa returned 
mirtars<-mieaa_all_gsea[mieaa_all_gsea$Category=='Target genes (miRTarBase)',]
# select columns mirnas, gene targets
all_targets<-mirtars[c('Subcategory', 'miRNAs.precursors')]


dim(all_targets)
colnames(all_targets)
library(data.table)

dt<-data.table(all_targets)
strsplit(all_targets$miRNAs.precursors, "\\; ")
all_targets_wide<-dcast(dt[, {x1 <- strsplit(miRNAs.precursors, "\\; "); c(list(unlist(x1)), 
          .SD[rep(seq_len(.N), lengths(x1))])}], Subcategory + miRNAs.precursors ~ V1, length)

all_targets_wide$miRNAs.precursors<-NULL
all_targets_long<-melt(all_targets_wide)
all_targets_long_true<-all_targets_long[all_targets_long$value==1, ]



### fix anticorrelation matrix
anticor_long<-melt(anticor)
anticor_long_ints<-anticor_long[anticor_long$value,]
colnames(anticor_long_ints)<-c('target_ensembl', 'mature_mirna_id', 'int')

symb<-get_symbols_vector(as.character(anticor_long_ints$target_ensembl))
anticor_long_ints$symbol<-symb
non_na<-(!is.na(symb)) & (!is.na(names(symb)))
## filter by the ones that returned syumbols 
anticor_long_ints_symbol<-anticor_long_ints[non_na,] 


anticor_long_ints
#symb_no_na<-symb_no_na[]
## WE LOST SOME GENES HERE 
length(symb);length(which(non_na))

anticor_long_ints_symbol$mature_mirna_id<-gsub('\\.', '-', anticor_long_ints_symbol$mature_mirna_id)

head(all_targets_long);
colnames(all_targets_long_true)<-c('symbol', 'mature_mirna_id', 'int')

head(anticor_long_ints_symbol)

### now filter all targets by what is actually anticorrelated 
merged_targets<-merge(all_targets_long_true, anticor_long_ints_symbol, by=c('mature_mirna_id', 'symbol'))
merged_targets

write.csv(merged_targets, paste0(mir_results_file, 'gene_targets_filtered_TSCAN.csv'))


ONT='BP'
#### Now run enrichment analysis only by the top results 
## rank/ order the list by the top mirnas --> a combination of pvalue and significance 
## i also get a pvalue for the interactions? 
enrich_go_mirnas=enrichGO(merged_targets$symbol, 
                       OrgDb='org.Hs.eg.db',
                       keyType = 'SYMBOL',
                       ont=ONT, 
                      pvalueCutoff = 0.05)
  

go_enrich<-enrich_go_mirnas
#View(enrich_go_mirnas@result)
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 30, 
        title = "GO Biological Pathways",
        font.size = 8)



##### plot results 
##
#DOSE::enrichResult
#  EnrichResult(mieaa_all_gsea)

  

df=mieaa_gsea_1
df$P.adjusted<-as.numeric(df$P.adjusted)
df$padj<-as.numeric(df$P.adjusted)


df$Observed<-as.numeric(df$Observed)

df_ord<-df[order(df$padj),]

df_ord<-df_ord[1:40,]
#df_ord$padj

mir_enrich_p<-ggplot(df_ord, aes(x=reorder(Subcategory, padj), y=Observed, fill=padj))+
  geom_bar(position='dodge', stat='identity')+
  coord_flip()



ggsave(paste0(mir_results_file, '_', Category, '_bar',  '.png'), height = 7, width=8)

mir_results_file_by_cat<-paste0(mir_results_file, '_', Category)

######### convert to enrichResult to use gsego functios


#install.packages("remotes")
#remotes::install_github("jmw86069/jamenrich")
library('multienrichjam')
library('clusterProfiler')

mieaa_gsea_1$P.adjusted<-as.numeric(mieaa_gsea_1$P.adjusted)


mieaa_gsea_1$keyColname=mieaa_gsea_1$Subcategory
mieaa_gsea_1_ord=mieaa_gsea_1[order(mieaa_gsea_1$P.adjusted),]
mieaa_gsea_1_ord$P.adjusted
#mieaa_gsea_1_ord_prob<-mieaa_gsea_1_ord
enr <- multienrichjam::enrichDF2enrichResult(as.data.frame(mieaa_gsea_1_ord),
                                             keyColname =  'Subcategory',
                                             geneColname ='miRNAs.precursors',
                                             pvalueColname = 'P.adjusted', 
                                             pvalueCutoff = 0.05)

# descriptionColname = "Subcategory",

enr
x2<-pairwise_termsim(enr)
N=15
p_emap<-emapplot(x2, showCategory =N )
p_emap

ggsave(paste0(mir_results_file_by_cat, '_conv_emap',  '.png'), height = 7, width=8)


p2<-dotplot(enr, showCategory=25)
p2
ggsave(paste0(mir_results_file_by_cat, '_conv_dotplot',  '.png'), height = 7, width=8)


### requires source('RNAseq enrichment.R') # TODO: MOVE TO A UTILS SCRIPT 
run_enrichment_plots(gse=enr, results_file=mir_results_file_by_cat)


############
###########




