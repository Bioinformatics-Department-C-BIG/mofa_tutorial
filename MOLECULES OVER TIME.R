


bl_f<-'ppmi/plots/p_BL_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_BL_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
bl<-read.csv(bl_f)

v08_f<-'ppmi/plots/p_V08_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
v08<-read.csv(v08_f)
bl<-bl[!is.na(bl$Description),]
tms<-list(bl=bl$Description, v08= v08$Description)
tms
out_com<-calculate.overlap()
out_com
venn.diagram(tms, 
             filename = paste0(outdir,prefix,'14_venn_diagramm.png'), output=FALSE)


v08_only<-v08[!(v08$Description %in% bl$Description),]
head(v08_only[order(v08_only$Least_value, decreasing = FALSE),]$Description, 50)
head(v08_only[order(v08_only$tot_rank, decreasing = FALSE),]$Description, 100)

v08_only
#### Markers over time:


ens_gene<-rownames(assay(se_filt_V08))[grep('ENSG00000115590', rownames(assay(se_filt_V08)))]
#ENSG00000196549
#ENSG00000151726
#ENSG00000138772
ens_gene<-rownames(assay(se_filt_V08))[grep('ENSG00000138772', rownames(assay(se_filt_V08)))]

se_filt_V08_pd<-se_filt_V08[,se_filt_V08$COHORT == 1]
se_filt_V08_pd$PD_MED_USE

med=se_filt_V08_pd$PD_MED_USE

se_filt_V08_pd<-se_filt_V08_pd[,se_filt_V08_pd$PD_MED_USE %in% c(0)]
se_filt_V08_pd$PD_MED_USE


se_filt_BL_pd<-se_filt_BL[,se_filt_BL$COHORT==1]
se_filt_BL$PATNO


se_filt_V08$PD_MED_USE
### AND ALSO changed STAGE 

df_v08<-assay(se_filt_V08_pd)

df_bl<-assay(se_filt_BL_pd)

v8_ens<-data.frame(df_v08[ens_gene,])
v8_ens$patno<-c(se_filt_V08_pd$PATNO)
bl_ens<-data.frame(df_bl[ens_gene,])
bl_ens$patno<-c(se_filt_BL_pd$PATNO)



merged_df<-merge(bl_ens, v8_ens, by='patno')


merged_df<-melt(merged_df )
ggpaired(merged_df, x = "variable", y = "value",
         color = "variable", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)




library(ggpubr)
# Create line plots of means
ggline(merged_df, x = 'variable', y = "value", 
       add = c("mean_sd", "jitter"))





