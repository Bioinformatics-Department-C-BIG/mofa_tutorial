#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
#install.packages('edgeR')
#BiocManager::install('limma')
#BiocManager::install('Glimma')

#### tutorial from: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html


library(limma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)

VISIT='V08'
source(paste0(script_dir, '/config.R' ))
source(paste0(script_dir,'/../bladder_cancer/preprocessing.R'))



##### START HERE WITH PROTEOMICS 
## TODO: SAVE AND LOAD 
# se_filt and vsn mat 

### THIS IS ALREADY filtered by cohort and VISIT 
# 
datalist<-loadRDS(prot_vsn_se_filt_file)
vsn_mat<-datalist[[1]]
se_filt<-datalist[[2]]

tmp<- assays(se_filt)[[1]]


if (run_vsn){
  protein_matrix<-vsn_mat
}else{
  protein_matrix<-tmp
  }

tmp
dim(proteomics)
COHORT<-se_filt$COHORT
dim(se_filt)
AGE=se_filt$AGE_SCALED
SEX=se_filt$SEX
COHORT

dim(se_filt)
formula_deseq2
AGE_AT_VISIT=AGE
formula_deseq_test<-'~AGE_AT_VISIT+SEX+COHORT'
design <- model.matrix(as.formula(formula_deseq_test) )
design <- model.matrix(~AGE_AT_VISIT+SEX+COHORT )
design <- model.matrix(~COHORT )
design <- model.matrix(~AGE_AT_VISIT+SEX+COHORT )


fit <- lmFit(protein_matrix, design = design)


#cont.matrix <- makeContrasts(B.NPS3vsNPS1=groupNPS1  - groupNPS3 ,levels=design)
#cont.matrix <- makeContrasts(pd_control=2-1,levels=design)

#fit.cont <- contrasts.fit(fit, cont.matrix)
#fit.cont <- eBayes(fit.cont, trend=TRUE)

fit.cont <- eBayes(fit, trend=TRUE)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)


library(ggplot2)
plotSA(fit)
nfeats<-dim(protein_matrix)[1]
nfeats
results_de<-topTable(fit.cont, number = nfeats, coef='COHORT' )
#results_de<-topTable(fit.cont, coef='COHORT' )
results_de
#FC= mean(condition_A_replicates)n/ mean(control_replicates)n   



  ### write out results 
# We want to highlight the significant genes. We can get this from decideTests.
#par(mfrow=c(1))
plotMD(fit.cont,coef=1,status=summa.fit[, 'COHORT' ], 
       values = c(-1, 1), hl.col=c("blue","red"))
dev.off()

fit.cont$ad
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
# TODO: create volcano plot with ggsave
#jpeg(paste0(outdir_s_p,'volcano.png'))
vp<-volcanoplot(fit.cont,coef=1,highlight=50,names=rownames(fit.cont$coefficients),
            main="pd_control"
            )
dev.off()

dir.create(outdir_s_p)

colnames(results_de)
ns_full<-table(se_filt$COHORT_DEFINITION)
ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])

# TODO: enhanced volcano
library(EnhancedVolcano)
ylim=max(-log10(results_de$adj.P.Val))+0.5
ylim
pvol<-EnhancedVolcano(results_de, 
                lab = rownames(results_de),
                x = 'logFC',
                y = 'adj.P.Val', 
                pCutoff = 10e-2,
                FCcutoff = 0.5, 
                ylim=c(0,ylim), 
                title='', 
                subtitle=ns
)
pvol
prefix='prot_'
fname<-paste0(outdir_s_p,'/EnhancedVolcano_edited_', prefix,'.jpeg')
ggsave(fname,pvol, width=6,height=8, dpi=300)
#ggsave(fname,pvol, width=6,height=8, dpi=300)


## Create a p-adjusted
## how many total proteins?
dim(vsn_mat)


dim(vsn_mat)[1]

#common_de<-intersect(all_sig_proteins,anova_results_oneway_significant)
#dir.create(outdir_s_p)
outdir_s_p_enrich<-paste0(outdir_s_p, '/enrichment/'); dir.create(outdir_s_p_enrich)
#write.csv(common_de, paste0(outdir_s_p, 'common_de.csv'))



#fit.cont_sig[common_de]
#




################### HEATMAPS  ############










#############################
#install.packages("OlinkAnalyze")
run_olink=FALSE
if (run_olink){
  
    
    library(OlinkAnalyze)
    #### ALSO TRY T-TEST
    data<-read_NPX(prot_files[1])
    data=prot_bl
    VISIT='BL'
    colnames(data)<- c( "PATNO" ,       "EVENT_ID"  ,   "Index"   ,     "OlinkID"    ,
                        "UniProt"     , "Assay"     ,   "MISSINGFREQ", "Panel"    , 
                        "PANEL_LOT_NR" ,"PLATEID"   ,   "QC_WARNING" ,  "LOD"    ,    
                        "NPX"   ,       "update_stamp")
    
    table(data$EVENT_ID)
    data_filt<-data[data$EVENT_ID %in% VISIT,]
    data_filt$EVENT_ID
    
    unique(data_filt$PATNO)
    outcome<-combined[,c('PATNO','COHORT', 'AGE_SCALED', 'SEX' )]
    outcome<-outcome[!duplicated(outcome$PATNO),]
    outcome<-outcome[outcome$COHORT %in% sel_coh,]
    
    ### now merge the specific data visit and metadata
    data_merged<-merge(data_filt, outcome, by='PATNO')
    data_merged$COHORT<-as.factor(data_merged$COHORT)
    
    data_merged$SampleID<-data_merged$PATNO
    ttest_results<-olink_ttest(df = data_merged,
                variable = 'COHORT')
    
    olink_wilcox_results<-olink_wilcox(df = data_merged,
                 variable = 'COHORT')
    
    olink_wilcox_results%>%
      dplyr::filter(Threshold == 'Significant')
    
    
    library(dplyr)
    data_merged$SEX=as.factor(data_merged$SEX)
    data_merged_rescaled<-data_merged
    data_merged_rescaled$NPX=data_merged_rescaled$NPX*1e6
    anova_results_oneway <- olink_anova(df = data_merged_rescaled, 
                                        variable = c('COHORT' ,'AGE_SCALED', 'SEX'))
    
    data_merged$COHORT=as.factor(data_merged$COHORT)
    anova_results_oneway <- olink_anova(df = data_merged, 
                                        variable = 'COHORT',
                                        covariates = c('AGE_SCALED', 'SEX'))
    
    anova_results_oneway_significant <- anova_results_oneway %>%
      dplyr::filter(Threshold == 'Significant')
    
    
    anova_results_oneway_significant_genes <- anova_results_oneway %>%
      dplyr::filter(Threshold == 'Significant') %>%
      dplyr::pull(Assay, Adjusted_pval)
    
    length(anova_results_oneway_significant)
    
    anova_results_oneway_significant[1:20]
    anova_results_oneway$Assay
    write.csv(anova_results_oneway,paste0(output_files, 'olink_de.csv' ))
    
    
    anova_posthoc_oneway_results <- olink_anova_posthoc(df = data_merged, 
                                        variable = 'COHORT',
                                        covariates = c('AGE_SCALED', 'SEX'))
    
    


}

T=0.05

################# ENRICHMENT - GSEA-GO #############
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'adj.P.Val'
order_statistic<-'signlog2pval'
order_statistic<-'logFC'
order_statistic<-'P.Value'
order_statistic<-'log2pval'
#order_statistic<-'pval_reverse'
#order_statistic<-'logFC'



#order_statistic<-'log2pval_not_adj' - NO RESULTS 
results_de$pval_reverse<- -results_de$P.Value

results_de$log2pval<- -log10(results_de$adj.P.Val) * results_de$logFC
results_de$abslog2pval<- abs(results_de$log2pval)

results_de$log2pval_not_adj<- -log10(results_de$P.Value) * results_de$logFC
results_de$signlog2pval<- -log10(results_de$adj.P.Val) * sign(results_de$logFC)


write.csv(results_de, paste0(outdir_s_p, 'results.csv'))


log2fol_T_overall<-0.1
padj_T_overall<-.05
results_de_signif<-mark_signficant(results_de,padj_T = padj_T_overall, log2fol_T = log2fol_T_overall, 
                            padj_name ='adj.P.Val',log2fc_name = 'logFC' , outdir_single = outdir_s_p  )

results_de_signif$abslog2pval

################### run gsea with anova ######################
gene_list1<-results_de[,order_statistic]
names(gene_list1)<-rownames(results_de)
gene_list_ord=gene_list1[order(-gene_list1)]
gene_list_ord
gene_list_limma_significant=rownames(results_de)[results_de$adj.P.Val<T]
gene_list_limma_significant_pval=rownames(results_de)[results_de$P.Value<T]


run_anova=FALSE
if (run_anova){
  order_statistic<-'statistic'
  gene_list2<-pull(anova_results_oneway, statistic)
  names(gene_list2)<-anova_results_oneway$Assay
  gene_list_ord=gene_list2[order(-gene_list2)]
  
  
}
length(gene_list_ord)


ONT='BP'
gene_list_ord


outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT, '_', order_statistic)
res_path<-paste0(outdir_s_p_enrich_file, 'gse.RDS')
if (file.exists(res_path)){
  gse_protein_full=loadRDS(res_path)
}else{
  
    
    #### TODO: ALSO RUN ENRICHMENT using ANOVA from 
    gse_protein_full <- clusterProfiler::gseGO(gene_list_ord, 
                                  ont=ONT, 
                                  keyType = 'SYMBOL', 
                                  OrgDb = 'org.Hs.eg.db', 
                                  pvalueCutoff  = pvalueCutoff)

  saveRDS(gse_protein_full, res_path)
}
dim(gse_protein_full@result)
hist(gse_protein_full@result$pvalue)
sig_ind<-gse_protein_full@result$pvalue<0.05
sig_ind
process_mirnas=TRUE

write.csv(as.data.frame(gse_protein_full@result), paste0(outdir_s_p_enrich_file, pvalueCutoff, '.csv'))

dim(gse_protein@result)
sig_gse_result<-gse_protein_full@result[gse_protein_full@result$pvalue<pvalueCutoff_sig,]
write.csv(as.data.frame(sig_gse_result), paste0(outdir_s_p_enrich_file,pvalueCutoff_sig ,'.csv'))

gse_protein_full[sig_ind,]
gse_protein<-gse_protein_full[sig_ind,]
enrich_plots<-run_enrichment_plots(gse=gse_protein_full,results_file=outdir_s_p_enrich_file , N_DOT=30, N_EMAP = 50)







### run enrichGo with anova

run_ORA=TRUE
#### IN ENRICHMENT WE filter by what is significant ONLY !! 
run_anova=FALSE
use_pval=FALSE
if (run_anova){
  order_statistic<-'statistic'
  anova_results_oneway$a
 anova_results_oneway_significant_en <- anova_results_oneway %>%
  dplyr::filter(Adjusted_pval<T)
  gene_list_anova_signficant<-pull(anova_results_oneway_significant_en, Assay)
  gene_list_ora=gene_list_anova_signficant
  print(length(gene_list_anova_signficant))
}else{
  gene_list_ora=gene_list_limma_significant
  if (use_pval){
    gene_list_ora=gene_list_limma_significant_pval
    
  }
  
  length(gene_list_limma_significant)
}

pvalueCutoff=1

gse_protein_enrich_full <- clusterProfiler::enrichGO(gene_list_ora, 
                                        ont=ONT, 
                                        keyType = 'SYMBOL', 
                                        OrgDb = 'org.Hs.eg.db', 
                                        pvalueCutoff  = pvalueCutoff)
## EMAP 15 is better for reporting - less crowded
process_mirnas=TRUE
#gse_protein_enrich_full<-gse_protein_enrich
outdir_s_p_enrich_file_ora<-paste0(outdir_s_p_enrich, ONT,  '_ora_T_', T,'_anova_', run_anova, 'pval_', use_pval )

gse_protein_enrich=write_filter_gse_results(gse_protein_enrich_full, outdir_s_p_enrich_file_ora, pvalueCutoff)
outdir_s_p_enrich_file_ora

#write.csv(as.data.frame(gse_protein_enrich@result), paste0(outdir_s_p_enrich_file_ora, '.csv'))
#write.csv(as.data.frame(gse_protein_enrich@result), paste0(outdir_s_p_enrich_file_ora, pvalueCutoff, '.csv'))
#pvalueCutoff_sig<-0.05
#sig_gse_result<-gse_protein_enrich@result[gse_protein_enrich@result$pvalue<pvalueCutoff_sig,]
#write.csv(as.data.frame(sig_gse_result), paste0(outdir_s_p_enrich_file_ora, pvalueCutoff_sig, '.csv'))
#dim(gse_protein_enrich@result)
#dim(sig_gse_result)


## does it contain nonsign pvalues? 
gse_protein_enrich@result[gse_protein_enrich@result$pvalue>0.05,]

enrich_plots<-run_enrichment_plots(gse=gse_protein_enrich,
                                   results_file=paste0(outdir_s_p_enrich_file_ora, 'enrich'), 
                                   N_DOT=30, N_EMAP = 15)

outdir_s_p_enrich_file_ora

dev.off()

gse_protein_enrich$p.adjust

### draw histogram 
df<-gse_protein_enrich@result
hist_p<-ggplot(df, aes(x=-log10(pvalue))) + 
  geom_histogram(color="black", fill="white")
hist_p

ggsave(paste0(outdir_s_p_enrich_file_ora, '_pval_hist.png'),plot=last_plot() )

int_params<-paste0(padj_paths, '_', VISIT, '_p_anova_',run_anova, 'pval_', use_pval )



