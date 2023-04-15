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
source('../bladder_cancer/preprocessing.R')
source(paste0(script_dir, '/config.R' ))


##### START HERE WITH PROTEOMICS 
tmp<- assays(se_filt)[[1]]
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

fit <- lmFit(vsn_mat, design = design)
#cont.matrix <- makeContrasts(B.NPS3vsNPS1=groupNPS1  - groupNPS3 ,levels=design)
#cont.matrix <- makeContrasts(pd_control=2-1,levels=design)

#fit.cont <- contrasts.fit(fit, cont.matrix)
#fit.cont <- eBayes(fit.cont, trend=TRUE)

fit.cont <- eBayes(fit, trend=TRUE)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)


library(ggplot2)
plotSA(fit)
nfeats<-dim(vsn_mat)[1]
results_de<-topTable(fit.cont, number = nfeats, coef='COHORT' )
#results_de<-topTable(fit.cont, coef='COHORT' )
results_de
#FC= mean(condition_A_replicates)n/ mean(control_replicates)n   



  ### write out results 
# We want to highlight the significant genes. We can get this from decideTests.
#par(mfrow=c(1))
plotMD(fit.cont,coef=1,status=summa.fit[,"pd_control"], 
       values = c(-1, 1), hl.col=c("blue","red"))
dev.off()

fit.cont$ad
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
# TODO: create volcano plot with ggsave
jpeg(paste0(outdir_s_p,'volcano.png'))
vp<-volcanoplot(fit.cont,coef=1,highlight=50,names=rownames(fit.cont$coefficients),
            main="pd_control"
            )
dev.off()
#ggsave(paste0(outdir_s_p,'volcano.png'), vp)

dir.create(outdir_s_p)


# TODO: enhanced volcano
library(EnhancedVolcano)
EnhancedVolcano(results_de)

## Create a p-adjusted
## how many total proteins?

fit.cont_sig<-fit.cont[fit.cont$p.value<0.05/dim(vsn_mat)[1],]
all_sig_proteins<-rownames(fit.cont_sig$p.value)

dim(vsn_mat)[1]

common_de<-intersect(all_sig_proteins,anova_results_oneway_significant)
#dir.create(outdir_s_p)
outdir_s_p_enrich<-paste0(outdir_s_p, '/enrichment/'); dir.create(outdir_s_p_enrich)
write.csv(common_de, paste0(outdir_s_p, 'common_de.csv'))

all_sig_proteins$

fit.cont_sig

#fit.cont_sig[common_de]
#











#############################
#install.packages("OlinkAnalyze")
library(OlinkAnalyze)
#### ALSO TRY T-TEST
data<-read_NPX(prot_files[1])
data=prot_bl

colnames(data)<- c( "PATNO" ,       "EVENT_ID"  ,   "Index"   ,     "OlinkID"    ,
                    "UniProt"     , "Assay"     ,   "MISSINGFREQ", "Panel"    , 
                    "PANEL_LOT_NR" ,"PLATEID"   ,   "QC_WARNING" ,  "LOD"    ,    
                    "NPX"   ,       "update_stamp")


outcome<-combined[,c('PATNO','COHORT', 'AGE_SCALED', 'SEX' )]
outcome<-outcome[!duplicated(outcome$PATNO),]
outcome<-outcome[outcome$COHORT %in% sel_coh,]
data_merged<-merge(data, outcome, by='PATNO')
data_merged$COHORT<-as.factor(data_merged$COHORT)

data_merged$SampleID<-data_merged$PATNO
ttest_results<-olink_ttest(df = data_merged,
            variable = 'COHORT')

ttest_results
library(dplyr)
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






T=0.05

################# ENRICHMENT - GSEA-GO #############
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'adj.P.Val'
order_statistic<-'signlog2pval'


#order_statistic<-'log2pval_not_adj' - NO RESULTS 


results_de$log2pval<- -log10(results_de$adj.P.Val) * results_de$logFC
results_de$log2pval_not_adj<- -log10(results_de$P.Value) * results_de$logFC
results_de$signlog2pval<- -log10(results_de$adj.P.Val) * sign(results_de$logFC)

################### run gsea with anova ######################
gene_list1<-results_de[,order_statistic]
names(gene_list1)<-rownames(results_de)
gene_list_ord=gene_list1[order(-gene_list)]

gene_list_limma_significant=rownames(results_de)[results_de$adj.P.Val<T]


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

### TODO: ALSO RUN ENRICHMENT using ANOVA from 
gse_protein <- clusterProfiler::gseGO(gene_list_ord, 
                              ont=ONT, 
                              keyType = 'SYMBOL', 
                              OrgDb = 'org.Hs.eg.db', 
                              pvalueCutoff  = 0.05)
gse@result

process_mirnas=TRUE
write.csv(as.data.frame(gse_protein@result), paste0(outdir_s_p_enrich_file, '.csv'))

outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT, '_', order_statistic)
enrich_plots<-run_enrichment_plots(gse=gse_protein,results_file=outdir_s_p_enrich_file , N_DOT=30, N_EMAP = 50)


### run enrichGo with anova

run_ORA=TRUE
#### IN ENRICHMENT WE filter by what is significant ONLY !! 
run_anova=FALSE
if (run_anova){
  order_statistic<-'statistic'
  anova_results_oneway$a
 anova_results_oneway_significant_en <- anova_results_oneway %>%
  dplyr::filter(Adjusted_pval<T)
  gene_list_anova_signficant<-pull(anova_results_oneway_significant_en, Assay)
  gene_list_ora=gene_list_anova_signficant
  length(gene_list_anova_signficant)
}else{
  gene_list_ora=gene_list_limma_significant
  length(gene_list_limma_significant)
}
gse_protein_enrich <- clusterProfiler::enrichGO(gene_list_ora, 
                                        ont=ONT, 
                                        keyType = 'SYMBOL', 
                                        OrgDb = 'org.Hs.eg.db', 
                                        pvalueCutoff  = 0.05)
## EMAP 15 is better for reporting - less crowded
outdir_s_p_enrich_file_ora<-paste0(outdir_s_p_enrich, ONT,  '_ora_T_', T,'_anova_', run_anova )
write.csv(as.data.frame(gse_protein_enrich@result), paste0(outdir_s_p_enrich_file_ora, '.csv'))
enrich_plots<-run_enrichment_plots(gse=gse_protein_enrich,
                                   results_file=paste0(outdir_s_p_enrich_file_ora, 'enrich'), 
                                   N_DOT=30, N_EMAP = 15)

outdir_s_p_enrich_file_ora








