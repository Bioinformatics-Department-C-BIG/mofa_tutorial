
############################# COMPARE PATHWAYS ###################################
##################################################################################
# for each visit across ALL modalities 
#### Load pathways from each of the modalities separately 


script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir,'ppmi/setup_os.R'))


#### Metascripts 
#('UpSetR')
library('UpSetR')
library('dplyr')

library('VennDiagram')
library(grid)
source(paste0(script_dir, 'ppmi/utils.R'))
source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))

process_mirnas=FALSE


VISIT='V08'; 
process_mirnas=TRUE; source(paste0(script_dir, 'ppmi/config.R' ))
outdir_mirs<-outdir_s; outdir_s
process_mirnas=FALSE; source(paste0(script_dir, 'ppmi/config.R' ))
outdir_rnas<-outdir_s; outdir_s


source(paste0(script_dir, 'ppmi/deseq_analysis_setup.R' ))
outdir_proteins<-outdir_s_p
padj_T_overall=0.05
log2fol_T_overall=0.1


## parameters for the enrichment- could be specified elsewhere?
run_anova=FALSE;use_pval=TRUE; 
run_ORA=FALSE; use_protein_pval=FALSE
padj_T=1;log2fol_T=0.00;order_by_metric<-'log2pval'; ONT='BP'



#outdir_s_p_enrich_file_ora<-paste0(outdir_proteins, '/enrichment/', 'BP_ora_T_0.05_anova_', run_anova,'pval_', use_pval)
outdir_s_p_enrich=paste0(outdir_proteins, '/enrichment/');order_statistic='logFC'
outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT,  '_', order_statistic, '_ora_', run_ORA,'ppval_', use_protein_pval, '_anova_', run_anova, 'pval_', use_pval )

results_file<-paste0(outdir_rnas, '/enrichment/', '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)



### obtain significant features from each one separately 
signif_rna<-read.csv(paste0(outdir_rnas,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))
signif_mirs<-read.csv(paste0(outdir_mirs,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))
signif_proteins<-read.csv(paste0(outdir_s_p,'/significant', padj_T_overall, '_',log2fol_T_overall, '.csv'))

### THE SIGNIFICANT GENES ARE TOO MANY.. filter them somehow..? 
signif_proteins



pvalueCutoff_sig=0.05
padj_paths<-0.05
pvalueCutoff=1; order_by_metric<-'log2pval'
results_file<-paste0(outdir_rnas,'/enrichment/', '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)
enrich_rnas_file<-paste0(results_file,pvalueCutoff ,'.csv')
#enrich_mirnas_file<-paste0(outdir_mirs,  '/enrichment/mirs_enrich__1_0_log2pval_GO Biological process (miRPathDB)',  '.csv')
enrich_mirnas_file<-paste0(outdir_mirs,  '/enrichment/GO Biological process (miRPathDB)/mirs_enrich__1_0_log2pvalGSEA600')
enrich_proteins_file<-paste0(outdir_s_p_enrich_file, pvalueCutoff, '.csv')


enrich_rna_single<-read.csv(enrich_rnas_file)
enrich_mirnas_single<-read.csv(paste0(enrich_mirnas_file,pvalueCutoff, '.csv'))
enrich_proteins_sigle<-read.csv(enrich_proteins_file)

enrich_rna_sig<-enrich_rna_single[enrich_rna_single$p.adjust<padj_paths,]; dim(enrich_rna_sig)[1]
enrich_mirnas_sig<-enrich_mirnas_single[enrich_mirnas_single$p.adjust<padj_paths,]; dim(enrich_mirnas_sig)[1]
enrich_proteins_sig<-enrich_proteins_sigle[enrich_proteins_sigle$p.adjust<padj_paths,]; dim(enrich_proteins_sig)[1]

common_paths<-intersect(enrich_rna_sig$Description,enrich_proteins_sig$Description )


## Create a list with all files 
listInput_all_mods_single<-list(rna=enrich_rna_sig$Description,
                                prot=enrich_proteins_sig$Description, 
                                mirnas=enrich_mirnas_sig$Description)
listInput<-listInput_all_mods_single




enrich_proteins_sig$p.adjust
get_ids(enrich_mirnas$ID[1])

length(unique(unlist(listInput_all_mods_single), use.names=FALSE))

res_overlap<-calculate.overlap(listInput)

intersection_all_three<-Reduce(intersect,listInput_all_mods)
int_params<-paste0(padj_paths, '_', VISIT, '_p_anova_',run_anova, 'pval_', use_pval )
write.csv(intersection_all_three, paste0(out_compare,'interesction_pathways' , int_params, '.csv') , row.names = FALSE)


library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


### UNION OF ALL single files 
venn.diagram(listInput,
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             cex=2.5,
             cat.cex=2.5,
             filename = paste0(out_compare,'all_modalities_', int_params ,'venn_diagramm.png'), output=TRUE)

############### COMBINE PVALUES #####

