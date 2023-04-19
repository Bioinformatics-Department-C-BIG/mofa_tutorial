

#### Metascripts 
#('UpSetR')
library('UpSetR')
library('dplyr')

library('VennDiagram')
library(grid)
source(paste0(script_dir, '/utils.R'))

process_mirnas=FALSE

### Table of samples from all visits 
#### FOR EACH MODALITY SEAPARETELY 
out_compare<-'ppmi/plots/single/compare/'

### this sets the outdirectory too
source(paste0(script_dir, '/config.R'))


dir.create(out_compare)

log2fol_T<-0.1;padj_T<-.05;

### TODO: filter by threshold here!! 
signif_file<-paste0('/significant', padj_T, '_',log2fol_T, '.csv')


#### Firstly print numbers of samples 
visits<-c('BL', 'V04', 'V06', 'V08')
all_vs_ps<-filter_se(se,visits, sel_coh )

meta<-colData(all_vs_ps)
get_stats<-meta[,c('PATNO', 'EVENT_ID')]
table_counts<-table(get_stats)
common_patients<-rownames(table_counts[rowSums(table_counts)==dim(table_counts)[2],])


list_all_vs<-split(get_stats, f = get_stats$EVENT_ID)



ns<-lapply(list_all_vs, function(x){
  length(unique(x$PATNO))
})
patient_lists<-lapply(list_all_vs, function(x){
  x$PATNO
})
upset(fromList(patient_lists))
counts_table<-table(get_stats)

## TODO: which are the common samples in all the lists


#View(c(meta$PATNO,meta$EVENT_ID))
library(plyr)
meta %>% 
  filter(EVENT_ID==visits[1])


visits=c('BL', 'V04', 'V06', 'V08')

outdir_s<-paste0(outdir_orig, '/single/', param_str_g_f, des)
outdir_s
if  (process_mirnas){
  outdir_all<-paste0(outdir_orig, '/single/', 'mirnas_',visits, '_',MIN_COUNT_M, '_coh_',sel_coh_s, '_',des)
  title_x='mirna'
}else{
  outdir_all<-paste0(outdir_orig, '/single/',  'rnas_',visits, '_',MIN_COUNT_G, '_coh_',sel_coh_s, '_', des)
  title_x='rna'
  
}


all_visits<-lapply(outdir_all, function(x) {
  
  outfile<-paste0(x, signif_file)
  if (file.exists(outfile)){
    return(as.data.frame(read.csv(outfile)))
  }else{
    print(paste0(outfile,' not available'))
  }
})


all_visits
padj_T<-0.01
log2fol_T<-0.1
list_of_mirs<-lapply(all_visits,function(df)
  if ('X' %in% colnames(df)){
    df$sign_lfc2 <- ifelse(df$padj <padj_T & abs(df$log2FoldChange) >log2fol_T , TRUE, FALSE);
    return(df[df$sign_lfc2,]$X)
  }
)
#list_of_mirs[[1]]
#merge(list_of_mirs, by='feats',all=TRUE)


listInput <- list(BL = list_of_mirs[[1]],
                  V04 =  list_of_mirs[[2]], 
                  V06 =  list_of_mirs[[3]],
                  V08 =  list_of_mirs[[4]])




jpeg(paste0(out_compare, title_x,'_upSet_diagram', padj_T, '_', log2fol_T, '.jpeg'), 
     res=200, width=800, height=500)
up<-upset(fromList(listInput), 
          sets.x.label = paste0(title_x," counts by Visit"))
up
dev.off()

#data_with_intersection <- listInput %>%
#  unite(col = "intersection", -c("entry"), sep = "")

calculate.overlap(listInput)

venn.diagram(listInput,   
             filename = paste0(out_compare,prefix,'14_venn_diagramm.png'), output=TRUE)




############################# COMPARE PATHWAYS ###################################
##################################################################################
# for each visit across ALL modalities 
#### Load pathways from each of the modalities separately 
VISIT='V08'; 
source(paste0(script_dir, '/config.R' ))
outdir_mirs<-paste0(outdir_orig, '/single/', 'mirnas_',VISIT, '_',MIN_COUNT_M, '_coh_',sel_coh_s, '_',des)
outdir_rnas<-paste0(outdir_orig, '/single/', 'rnas_',VISIT, '_',MIN_COUNT_G, '_coh_',sel_coh_s, '_',des)
#outdir_proteins<-paste0(outdir_orig, '/single/proteomics_', VISIT,'_coh_', sel_coh_s, '_', des, '/' )
outdir_proteins<-outdir_s_p



## parameters for the enrichment- could be specified elsewhere?
run_anova=FALSE;use_pval=TRUE; 
padj_T=1;log2fol_T=0.00;order_by_metric<-'log2pval'; ONT='BP'

outdir_s_p_enrich_file_ora<-paste0(outdir_proteins, '/enrichment/', 'BP_ora_T_0.05_anova_', run_anova,'pval_', use_pval)
results_file<-paste0(outdir_rnas, '/enrichment/', '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)



### obtain significant features from each one separately 
log2fol_T_overall=0.1
signif_rna<-read.csv(paste0(outdir_rnas,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))
signif_mirs<-read.csv(paste0(outdir_mirs,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))
signif_proteins<-read.csv(paste0(outdir_s_p,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))

### THE SIGNIFICANT GENES ARE TOO MANY.. filter them somehow..? 






padj_paths<-0.05
enrich_rnas_file<-paste0(results_file, '.csv')
enrich_mirnas_file<-paste0(outdir_mirs,  '/enrichment/mirs_enrich__1_0_log2pval_GO Biological process (miRPathDB)',  '.csv')
enrich_mirnas_file<-paste0(outdir_mirs,  '/enrichment/GO Biological process (miRPathDB)/mirs_enrich__1_0_log2pval',  '.csv')

enrich_proteins_file<-paste0(outdir_s_p_enrich_file_ora,'.csv')
enrich_rna<-read.csv(enrich_rnas_file)
enrich_mirnas<-read.csv(enrich_mirnas_file)
enrich_proteins<-read.csv(enrich_proteins_file)



enrich_rna_sig<-enrich_rna[enrich_rna$p.adjust<padj_paths,]
enrich_mirnas_sig<-enrich_mirnas[enrich_mirnas$p.adjust<padj_paths,]
dim(enrich_mirnas_sig)
enrich_proteins_sig<-enrich_proteins[enrich_proteins$p.adjust<padj_paths,]

common_paths<-intersect(enrich_rna_sig$Description,enrich_proteins_sig$Description )


listInput_all_mods<-list(rna=enrich_rna_sig$Description,
                         prot=enrich_proteins_sig$Description, 
                         mirnas=enrich_mirnas_sig$Description)



listInput<-listInput_all_mods
res_overlap<-calculate.overlap(listInput)

intersection_all_three<-Reduce(intersect,listInput_all_mods)
int_params<-paste0(padj_paths, '_', VISIT, '_p_anova_',run_anova, 'pval_', use_pval )
write.csv(intersection_all_three, paste0(out_compare,'interesction_pathways' , int_params, '.csv') , row.names = FALSE)
unique_rna<-enrich_rna %>% 
  filter(Description %in% res_overlap$a1)
dim(unique_rna)
lapply(res_overlap, length)

write.csv(unique_rna,
          paste0(out_compare,'unique_rna_' , int_params, '.csv') , row.names = FALSE)

write.csv(enrich_proteins %>% 
            filter(Description %in% res_overlap$a3), 
          paste0(out_compare,'unique_prot_' , int_params, '.csv') , row.names = FALSE)


write.csv(enrich_mirnas %>% 
            filter(Description %in%res_overlap$a7),
          paste0(out_compare,'unique_mirs_' , int_params, '.csv') , row.names = FALSE)



library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


venn.diagram(listInput,
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             cex=2.5,
              
             
             
             filename = paste0(out_compare,'all_modalities_', int_params ,'venn_diagramm.png'), output=TRUE)


