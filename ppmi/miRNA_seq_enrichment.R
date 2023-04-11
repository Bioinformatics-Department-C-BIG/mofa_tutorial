
#install.packages('VennDiagram')
library(rbioapi)
library('VennDiagram')
library(enrichplot)

library('GOfuncR')
require(DOSE)
library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library('org.Hs.eg.db')
#install.packages('ggridges')
require(ggridges)


source('ppmi/deseq_analysis.R')

#order_by_metric<-'abslog2pval'
#order_by_metric<-'abslog2pval'
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




######################################################################
### ANOTHER WAY to obtain targets 
##
##

#BiocManager::install("targetscan.Hs.eg.db")
### first retrieve all targets that mieaa returned 
####################

##### wgere to draw targets from????? 
mir_results_file_anticor<-paste0(mir_results_file, '_anticor_')



mirtars<-mieaa_all_gsea[mieaa_all_gsea$Category=='Target genes (miRTarBase)',]
# select columns mirnas, gene targets
all_targets<-mirtars[c('Subcategory', 'miRNAs.precursors')]

library(data.table)
dt<-data.table(all_targets)
strsplit(all_targets$miRNAs.precursors, "\\; ")
all_targets_wide<-dcast(dt[, {x1 <- strsplit(miRNAs.precursors, "\\; "); c(list(unlist(x1)), 
          .SD[rep(seq_len(.N), lengths(x1))])}], Subcategory + miRNAs.precursors ~ V1, length)

all_targets_wide$miRNAs.precursors<-NULL
all_targets_long<-melt(all_targets_wide)
all_targets_long_true<-all_targets_long[all_targets_long$value==1, ]





head(all_targets_long);
colnames(all_targets_long_true)<-c('symbol', 'mature_mirna_id', 'int')


##################### METHOD 2: RANK BY ANTI COR ######################
#################


########### fix anticorrelation matrix
cor_results_long_all<-melt(cor_results, varnames = c('target_ensembl', 'mature_mirna_id'), value.name = 'cor' )
#cor_results_long_all<-cor_results_long
# could filter here
### filter using binary threhsold 

T_cor=-0.3
T_cor=0

#T_cor=1

mir_results_file_anticor=paste0(mir_results_file, '_anticor_T_cor_', T_cor)

cor_results_long_ints<-cor_results_long_all[cor_results_long_all$cor<=T_cor,]


#cor_results_long_ints<-cor_results_long_all[cor_results_long_all$cor<=-0.3,]
cor_results_long<-cor_results_long_ints


symb<-get_symbols_vector(as.character(cor_results_long$target_ensembl))
cor_results_long$symbol<-symb
non_na<-(!is.na(symb)) & (!is.na(names(symb)))
## filter by the ones that returned syumbols 
cor_results_long_symbol<-cor_results_long[non_na,] 
head(cor_results_long_symbol)
## WE LOST SOME GENES HERE 
length(symb);length(which(non_na))
cor_results_long_symbol$mature_mirna_id<-gsub('\\.', '-', cor_results_long_symbol$mature_mirna_id)

hist(cor_results_long_symbol$cor)

##### ALL POSSIBLE TARGETS ################################

colnames(all_targets_long_true)<-c('symbol', 'mature_mirna_id', 'int')



##################### merge all possible targets with correlation values 
###
merged_targets<-merge(all_targets_long_true, cor_results_long_symbol, by=c('mature_mirna_id', 'symbol'))
#hist(merged_targets$cor)



gene_list_metric<-DataFrame(order_by_metric=gene_list)
gene_list_metric$mature_mirna_id=rownames(gene_list_metric)
colnames(gene_list_metric)


merged_targets_metric<-merge(merged_targets,gene_list_metric,  by=c('mature_mirna_id'))
dim(merged_targets_metric)
dim(merged_targets)

order_metric='log2pval_negcor'
merged_targets_metric[,order_metric]<-merged_targets_metric$cor * -1 * merged_targets_metric$order_by_metric

#merged_targets_metric[,order_metric]<-merged_targets_metric$cor * -1 * merged_targets_metric$order_by_metric


hist(merged_targets_metric[,order_metric])


write.csv(merged_targets_metric, paste0(mir_results_file_anticor, 'gene_targets_filtered.csv'))


gene_list_targets<-merged_targets_metric[,order_metric]
names(gene_list_targets)<-merged_targets_metric$target_ensembl
gene_list_targets_ord<-gene_list_targets[order(-gene_list_targets)]


ONT='BP'
#### Now run enrichment analysis only by the top results 
## rank/ order the list by the top mirnas --> a combination of pvalue and significance 
## i also get a pvalue for the interactions? 
gse_mirnas= clusterProfiler::gseGO(gene_list_targets_ord, 
                                                 ont=ONT, 
                                                 keyType = 'ENSEMBL', 
                                                 OrgDb = 'org.Hs.eg.db', 
                                                 pvalueCutoff  = 0.05)




  
  
  
  

#View(enrich_go_mirnas@result)
write.csv(gse_mirnas@result, paste0(mir_results_file_anticor, 'results.csv'))

#ggsave(paste0(mir_results_file_anticor, '_',T_cor, '_barplot',  '.jpeg'), width=8, height=7)

#run_enrichment_plots(gse_mirnas, )
######################### RANKED BY NEGATIVE-CORELATION #########################



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

######### convert to enrichResult to use gsego functions ##########################
###################################################################################


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




