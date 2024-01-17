


## compare cell correction 

## corr cell is corrected by neutrophil score 
# 1. how many new paths/lost paths with correction
# 2. which paths have remained after corrections? 
# 3. what is the decrease in pvalue for the changed paths?

de_file_corr1 = paste0(outdir_s_corr1, '/results_df.csv')
results_file_enrich_corr1=paste0(outdir_s_corr1, '/enr/', prefix, enrich_params, '1.csv')

de_file_corr0 = paste0(outdir_s_corr0, '/results_df.csv')
results_file_enrich_corr0=paste0(outdir_s_corr0, '/enr/', prefix, enrich_params, '1.csv')


de_file_corr1_DF<-read.csv(de_file_corr1, row.names=1) 
de_file_corr0_DF<-read.csv(de_file_corr0, row.names=1) 

results_enrich_corr1_DF<-read.csv(results_file_enrich_corr1)
results_enrich_corr0_DF<-read.csv(results_file_enrich_corr0)


results_enrich_corr0_sig<-results_enrich_corr0 %>% dplyr::filter(p.adjust<0.05)

dim(results_enrich_corr0 %>% dplyr::filter(p.adjust<0.05))
results_enrich_corr1_sig<-results_enrich_corr1 %>% dplyr::filter(p.adjust<0.05)


intersect(results_enrich_corr1_sig$Description, results_enrich_corr0_sig$Description)
deseq_params_all = outdir_s
deseq_all_list

gl1
# get ordered gene list from de results 
# 
gl1<-get_ordered_gene_list(de_file_corr1_DF,  order_by_metric, padj_T=1, log2fol_T=0 ); # for corrected
names(gl1)<-gsub('\\..*', '',names(gl1)) # rename genelist  ens names

gl0<-get_ordered_gene_list(de_file_corr0_DF,  order_by_metric, padj_T=1, log2fol_T=0 );# and uncorected corrected
names(gl0)<-gsub('\\..*', '',names(gl0))# rename genelist ens names

names(gl1)
dir.create(paste0(deseq_params_all, '/enr/'))
enrich_compare_path=paste0(deseq_params_all, '/enr/','compare')


##
gse_compare<-compareCluster(geneClusters = list(corr=gl1, no_corr=gl0 ), 
                            fun = "gseGO", 
                            OrgDb='org.Hs.eg.db', 
                            ont=ONT, 
                            keyType = 'ENSEMBL') 


plot_enrich_compare(gse_compare,paste0(enrich_compare_path), N_EMAP = 30)


