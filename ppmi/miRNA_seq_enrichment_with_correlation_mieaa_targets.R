

#### Another enrichment analysis using correlation to rank genes
### ANOTHER WAY to obtain targets 
##
##
library(data.table)
#BiocManager::install('targetscan.Hs.eg.db')
library(targetscan.Hs.eg.db)
library(org.Hs.eg.db)

use_anticor=TRUE


log2fol_T=0.0

#BiocManager::install("targetscan.Hs.eg.db")
################## MIEAA TARGETS 
### first retrieve all targets that mieaa returned 
## Return a 2d table with mirnas in one column and possible gene targets in the other
mirtars<-mieaa_all_gsea[mieaa_all_gsea$Category=='Target genes (miRTarBase)',]
# select columns mirnas, gene targets
all_targets<-mirtars[c('Subcategory', 'miRNAs.precursors')]

dim(all_targets)
colnames(all_targets)

dt<-data.table(all_targets)
strsplit(all_targets$miRNAs.precursors, "\\; ")
all_targets_wide<-dcast(dt[, {x1 <- strsplit(miRNAs.precursors, "\\; "); c(list(unlist(x1)), 
                               .SD[rep(seq_len(.N), lengths(x1))])}], Subcategory + miRNAs.precursors ~ V1, length)
all_targets_wide$miRNAs.precursors<-NULL
all_targets_long<-melt(all_targets_wide)
all_targets_long_true<-all_targets_long[all_targets_long$value==1, ]
colnames(all_targets_long_true)


colnames(all_targets_long_true)<-c('symbol', 'mature_mirna_id', 'int')

### method1: returns all_targets_long_true: all targets from mieaa 



##############
#### Targets 2: TARGET SCAN ############
run_tscan=FALSE
if (run_tscan){
  

        all_target_scan<-mget(names(mirs), revmap(targetscan.Hs.egTARGETS))
        tars<-all_target_scan[[1]]
        
        length(all_target_scan)
        
        named_list<-map2(all_target_scan, names(all_target_scan),~cbind(.x, ID=.y)) 
        
        
        names1<-head(unlist(all_target_scan3))
        nlist_stack<-as.data.frame(do.call(rbind,named_list))
        nlist_stack
        colnames(nlist_stack)<-c('entrezid', 'mature_mirna_id')
        mirs_nohs<- gsub('hsa-', '', mirs)
        mirs_nohs
        nlist_stack_filt<-nlist_stack[nlist_stack$mature_mirna_id %in% mirs_nohs,]
        ids_ens<-mapIds(org.Hs.eg.db, keys = nlist_stack_filt$entrezid, keytype="ENTREZID", column = "ENSEMBL")
        ids_ens_df<-DataFrame(unlist(ids_ens)); 
        ids_ens_df$entrezid=rownames(ids_ens_df)
        colnames(ids_ens_df)=c('target_ensembl', 'entrezid')
        
        merged_gene_tars<-merge(ids_ens_df@listData, DataFrame(nlist_stack_filt), by='entrezid' )
        merged_gene_tars_unique<-unique(merged_gene_tars)
        colnames(merged_gene_tars_unique)
        }


########### fix anticorrelation matrix #####################
### Filter by anticorrelated?? 

T_cor=-0.1 ## if we dont then everything comes out at the same value... 
TOP_GN_all=0.9; TOP_MN_all=0.9
cor_results_read<-read.csv(paste0(outdir_s, '/',   TOP_GN_all, '_', TOP_MN_all, '_','cor_results.csv')) ## not sure try it 
cor_results_long_all<-melt(cor_results_read, varnames = c('target_ensembl', 'mature_mirna_id'), value.name = 'cor' )
colnames(cor_results_long_all)
colnames(cor_results_long_all)<-c('target_ensembl', 'mature_mirna_id', 'cor')
#cor_results_long_all<-cor_results_long
# could filter here
### filter using binary threhsold 

cor_results_long_ints<-cor_results_long_all[cor_results_long_all$cor<=T_cor,]

cor_results_long_ints
#cor_results_long_ints<-cor_results_long_all[cor_results_long_all$cor<=-0.3,]
cor_results_long<-cor_results_long_ints
turn_to_symb=TRUE
if (turn_to_symb){
  symb<-get_symbols_vector(as.character(cor_results_long$target_ensembl))
  cor_results_long$symbol<-symb
  non_na<-(!is.na(symb)) & (!is.na(names(symb)))
  ## filter by the ones that returned syumbols 
  cor_results_long_symbol<-cor_results_long[non_na,] 
  head(cor_results_long_symbol)
  length(symb);length(which(non_na))
  cor_results_long_symbol$mature_mirna_id<-gsub('\\.', '-', cor_results_long_symbol$mature_mirna_id)
  hist(cor_results_long_symbol$cor)
  
  
}

hist(cor_results_long_symbol$cor)
hist(cor_results_long$cor)

##################### merge all possible targets with correlation values 
###
head(cor_results_long_symbol)
head(all_targets_long_true)

merged_targets<-merge(all_targets_long_true, cor_results_long_symbol, by=c('mature_mirna_id', 'symbol'))


if (run_tscan){
      merged_gene_tars_unique$mature_mirna_id<-gsub('-', '.', merged_gene_tars_unique$mature_mirna_id)
      cor_results_long$mature_mirna_id<- gsub('hsa.', '', cor_results_long$mature_mirna_id)
      head(merged_gene_tars_unique); head(cor_results_long)
      
      
      merged_targets_scan<-merge(merged_gene_tars_unique, cor_results_long, by=c('mature_mirna_id', 'target_ensembl'))
      
      merged_targets<-merged_targets_scan 
      
      
      names(gene_list)<-gsub('-', '.', names(gene_list))
      names(gene_list)<- gsub('hsa.', '', names(gene_list))
      
}




#### Fix the gene list to match the targets list 

gene_list_metric<-DataFrame(order_by_metric=gene_list)
gene_list_metric$mature_mirna_id=rownames(gene_list_metric)
colnames(gene_list_metric)


merged_targets_metric<-merge(merged_targets,gene_list_metric,  by=c('mature_mirna_id'))
dim(merged_targets_metric)
dim(merged_targets)


order_metric='log2pval_negcor'
merged_targets_metric[,order_metric]<- merged_targets_metric$cor * -1 * merged_targets_metric$order_by_metric
merged_targets_metric

#merged_targets_metric[,order_metric]<-merged_targets_metric$cor * -1 * merged_targets_metric$order_by_metric
############# REMOVE DUPLICATED ###################
## manual check to see that 
## Here we should expect high negative correlations to have positive metric (well unless fc is negative !! )

mer_tars_ord<-merged_targets_metric[order(-merged_targets_metric$log2pval_negcor),]
remove_dup_g=FALSE

if (remove_dup_g){
  mer_tars_ord_no_dup<-mer_tars_ord[!duplicated(mer_tars_ord$target_ensembl),]
  
}else{
  mer_tars_ord_no_dup<-mer_tars_ord
  
}

### 

hist(mer_tars_ord_no_dup[,order_metric])
dir.create(paste0(outdir_enrich, '/anticor/'))
mir_results_file_anticor=paste0(outdir_enrich, '/T_cor_',enrich_params  , T_cor, 'dup', remove_dup_g )


write.csv(mer_tars_ord_no_dup, paste0(mir_results_file_anticor, 'gene_targets_filtered_',run_tscan, 'LAT', '.csv'))



##### ORDER THE GENE LIST ############################
gene_list_targets<-mer_tars_ord_no_dup[,order_metric]
names(gene_list_targets)<-mer_tars_ord_no_dup$target_ensembl
gene_list_targets_ord<-gene_list_targets[order(-gene_list_targets)]
write.csv(merged_targets, paste0(mir_results_file_anticor, 'gene_targets_filtered_TSCAN.csv'))


ONT='BP'

gene_list_targets_ord
### now we are filtering and using a targeted set of genes so run gse go!! 
gse_mirnas <- clusterProfiler::gseGO(gene_list_targets_ord, 
                                   ont=ONT, 
                                   keyType = 'ENSEMBL', 
                                   OrgDb = 'org.Hs.eg.db', 
                                   pvalueCutoff  = 0.05)






#View(enrich_go_mirnas@result)
write.csv(gse_mirnas@result, paste0(mir_results_file_anticor, 'results.csv'))

#ggsave(paste0(mir_results_file_anticor, '_',T_cor, '_barplot',  '.jpeg'), width=8, height=7)

run_enrichment_plots(gse_mirnas, mir_results_file_anticor)
######################### RANKED BY NEGATIVE-CORELATION #########################

dotplot(gse_mirnas, showCategory=30) + ggtitle("dotplot for ORA")


gse_mirnas@result$p.adjust


