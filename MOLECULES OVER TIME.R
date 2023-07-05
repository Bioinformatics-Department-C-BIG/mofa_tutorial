
library(org.Hs.eg.db)


## Load baseline and v08 ####
# run deseqvst
bl_f<-'ppmi/plots/p_BL_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_BL_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
bl<-read.csv(bl_f)

v08_f<-'ppmi/plots/p_V08_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
v08<-read.csv(v08_f)
bl<-bl[!is.na(bl$Description),]
bl_sig<-bl[bl$Least_value<0.01,]

unique(bl_sig$Description)
tms<-list(bl=bl$Description, v08= v08$Description)

venn.diagram(tms, 
             filename = paste0(outdir,prefix,'14_venn_diagramm.png'), output=FALSE)


v08_only<-v08[!(v08$Description %in% bl$Description),]
head(v08_only[order(v08_only$Least_value, decreasing = FALSE),]$Description, 50)
head(v08_only[order(v08_only$tot_rank, decreasing = FALSE),]$Description, 100)

v08_only
#### Markers over time:
#### 1. Obtain the markers here 
fn_sel=2
view='miRNA'; process_mirnas=TRUE
source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))

ws<-get_weights(MOFAobject, views = view, factors=sel_factors[fn_sel])[[1]]


ws
ws_high<-ws[abs(ws)>quantile(ws, 0.95) & abs(ws)<quantile(ws, 0.96) ,]
ws_high<-ws[abs(ws)>quantile(ws, 0.85) & abs(ws)<quantile(ws, 0.86) ,]
ws_high<-ws[abs(ws)>quantile(ws, 0.95) & abs(ws)<quantile(ws, 0.96) ,]
ws_high<-ws[abs(ws)>quantile(ws, 0.999999999),]
ws_high<-ws[abs(ws)>quantile(ws, 0.99999996),][1:12]

#ws_high<-ws_high[1:15]
ws_high
#ws_high<-ws_high[abs(ws_high)<quantile(ws_high, 0.9)]






if (view=='RNA'){
  ens_genes<-rownames(assay(se_filt_V08))[grep(paste0(names(ws_high), collapse='|'), rownames(assay(se_filt_V08)))]
  feat_names=ens_genes
  
}else{
  feat_names=names(ws_high)
  
}
se_filt_V08_pd<-se_filt_V08[,se_filt_V08$COHORT == 1]
se_filt_V08_pd<-se_filt_V08_pd[,se_filt_V08_pd$PATNO %in% common]
dim(se_filt_V08_pd)
dim(se_filt_V08_pd)

se_filt_V08_pd$PD_MED_USE
med=se_filt_V08_pd$PD_MED_USE
se_filt_BL_pd<-se_filt_BL[,se_filt_BL$COHORT==1]


### AND ALSO changed STAGE 

df_v08<-assay(se_filt_V08_pd)
df_bl<-assay(se_filt_BL_pd)
df_v08<-cpm(df_v08, normalized.lib.sizes=TRUE)
df_bl<-cpm(df_bl, normalized.lib.sizes=TRUE)



## outliers

clip_outliers<-function(df1){
  #'
  #' @param 
  #'
  #'
  df1.quantiles <- apply(df1, 1, function(x, prob=0.99) { quantile(x, prob, names=F) })
  for (i in 1:dim(df1)[1]){
    df1[i,][ df1[i,]> df1.quantiles[i] ]<- df1.quantiles[i]
  }
  
  return(df1)
}

df_bl<- clip_outliers(df_bl)
df_v08<- clip_outliers(df_v08)






######### PLOT molecular markers 


unique(colnames(df_v08))
colnames(df_v08)
rownames(df_v08)
rownames(df_v08) %in% ens_genes
dim(df_v08)

v8_ens<-t(df_v08[rownames(df_v08) %in% feat_names,])
v8_ens=data.frame(v8_ens,patno=c(se_filt_V08_pd$PATNO), PATNO_EVENT_ID=c(se_filt_V08_pd$PATNO_EVENT_ID))
bl_ens<-t(df_bl[rownames(df_bl) %in% feat_names,])
bl_ens=data.frame(bl_ens,patno=c(se_filt_BL_pd$PATNO), PATNO_EVENT_ID=c(se_filt_BL_pd$PATNO_EVENT_ID))
se_filt_BL_pd$PATNO_EVENT_ID


v8_ens
### MELT and MERGE 
v8_melt<-melt(v8_ens)
v8_melt$VISIT<-'V08'

bl_melt<-melt(bl_ens); bl_melt$VISIT='BL'
merged_melt<-rbind(bl_melt, v8_melt)
merged_melt
#merged_df<-merge(bl_ens, v8_ens, by='patno')
#merged_df

ens<-gsub('\\..*', '',merged_melt$variable)




#
### ### NOW match factors to samples
Z <- get_factors(MOFAobject)[[1]]
Z1<-Z[,sel_factors[fn_sel]]
patnos_z1<-gsub('\\_.*', '', names(Z1))
Z1_matched<-Z1[match(merged_melt$patno,patnos_z1) ]

png(paste0(outdir, '/trajectories/','hist_', sel_factors[fn_sel],'_', view,  '.jpeg'))
hist(Z1_matched);
dev.off()

quantile(Z1_matched,0.85, na.rm=TRUE)
Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,0.2, na.rm=TRUE))
Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,0.85, na.rm=TRUE))
Z1_grouping
dim(merged_melt)
merged_melt$grouping<-Z1_grouping


if (view=='RNA'){
  symb<-get_symbols_vector(ens)
  merged_melt$symbol<-symb
}else{
  merged_melt$symbol<-merged_melt$variable
  
}
merged_melt$variable
merged_melt
na_ps<-unique(merged_melt[!is.na(merged_melt$grouping),]$patno)
na_ps
merged_melt_filt<-merged_melt[merged_melt$patno %in% na_ps, ]
merged_melt_filt

### EDIT GROUP TO BE PATNO and GROUP!!
# OR JUST GROUP  
merged_melt_filt$VISIT








#### 
ggplot(data = merged_melt_filt, aes(x = grouping, y = value)) + 
  #geom_point(aes(col=factor(VISIT)), size = 2) +
  #geom_line(aes(group=patno, col=grouping)) +
  geom_boxplot(aes(fill=VISIT ))+
  geom_line(aes(group=patno ))+
  
  #geom_line(aes(group=patno, col=VISIT)) +
  
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  facet_wrap(. ~ symbol, scales='free_y') +
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw() + 
  
  
  stat_compare_means(comparisons = list(c("BL", "V08")), 
                     label = "p.format", method = "wilcox.test", tip.length = 0)

warnings()
ggsave(paste0(outdir, '/trajectories/', sel_factors[fn_sel],'_', view,  '.jpeg'), 
       width=12, height=12)

merged_melt_filt$VISIT
merged_melt_filt$VISIT<-as.numeric(factor(merged_melt_filt$VISIT))
















stat.test2 <- merged_melt_filt %>%
  group_by(grouping) %>%
  t_test(value ~ VISIT, p.adjust.method = "bonferroni")
stat.test2

##############

### Plot clinical variables 

na_ps
### AND ALSO changed STAGE 

se_V08_cl<-se_filt_V08_pd[,se_filt_V08_pd$PATNO%in% na_ps ]
se_bl_cl<-se_filt_BL_pd[,se_filt_BL_pd$PATNO%in% na_ps ]


se_V08_cl$NP4TOT
sel_factors[fn_sel]
to_sel_cor<-names(cors_all[names(sel_factors)[fn_sel],][cors_all[names(sel_factors)[fn_sel],]>1.7])
to_sel<- c('PATNO', to_sel_cor)
to_sel<-intersect(names(colData(se_V08_cl)), to_sel)
df_V08_cl<-colData(se_V08_cl)[,to_sel]
df_bl_cl<-colData(se_bl_cl)[,to_sel]

df_bl_cl$MCAVFNUM

df_V08_cl_ml<-reshape2::melt(df_V08_cl,value.name = to_sel); df_V08_cl_ml$VISIT='V08'
df_bl_cl_ml<-reshape2::melt(df_bl_cl, value.name=to_sel); df_bl_cl_ml$VISIT='BL'

df_V08_cl_ml


## Which variables are corelated? 


merged_melt_cl<-rbind(df_bl_cl_ml,df_V08_cl_ml)

patnos_z1<-gsub('\\_.*', '', names(Z1))
patnos_z1
merged_melt_cl$PATNO
Z1_matched<-Z1[match(merged_melt_cl$PATNO,patnos_z1) ]
Z1
hist(Z1_matched);
quantile(Z1_matched,0.9, na.rm=TRUE)
Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,0.9, na.rm=TRUE))
Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,0.2, na.rm=TRUE))
Z1_grouping
dim(merged_melt)
merged_melt_cl$grouping<-Z1_grouping

to_sel

### for categorical 
table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'],merged_melt_cl$grouping )
table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'], )
merged_melt_cl$NP3SPCH
ggplot(data = merged_melt_cl, aes( x=factor(grouping), 
                                   fill = factor(PDSTATE) )) + 
  geom_bar()+
  facet_wrap(. ~ VISIT, scales='free_y') 
  
#theme_bw() 
to_sel

### for continous 
is.numeric(merged_melt_cl$LAST_UPDATE_M1)
cov_to_plot %in% to_sel

table(merged_melt_cl$PDSTATE)
merged_melt_cl3$PDSTATE
merged_melt_cl3<-merged_melt_cl[merged_melt_cl$PDSTATE %in% c('OFF', ''),]
merged_melt_cl2<-merged_melt_cl3

for (cov_to_plot in to_sel){
  
  
    if (is.numeric(merged_melt_cl3[, cov_to_plot])){
      
    
    ggplot(data = merged_melt_cl3, aes_string(x = 'grouping', y =cov_to_plot  )) + 
      geom_boxplot(aes(fill=VISIT ))+
    #facet_wrap(. ~ symbol, scales='free_y') +
        theme_bw() 
      
      ggsave(paste0(outdir, '/trajectories/', 'cl_var', cov_to_plot, sel_factors[fn_sel],to_sel[2]  , '.jpeg'))
       merged_melt_cl2[,cov_to_plot]<-factor(merged_melt_cl2[,cov_to_plot])
       
      ggplot(data = merged_melt_cl2, aes_string( x='grouping', 
                                         fill = cov_to_plot) ) + 
        geom_bar()+
        facet_wrap(. ~ VISIT, scales='free_y') 
      
      ggsave(paste0(outdir, '/trajectories/', 'cl_var_discrete', cov_to_plot, sel_factors[fn_sel],to_sel[2]  , '.jpeg'))
      
    }
}

ggsave(paste0(outdir, '/trajectories/', sel_factors[fn_sel],to_sel[2]  , '.jpeg'))




table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3TOT'],merged_melt_cl$grouping )


ggplot(data = merged_melt_cl, aes_string(x = 'VISIT', y ='NP3SPCH'  )) + 
  geom_point(aes(col=factor(VISIT)), size = 2) +
  geom_line(aes(group=PATNO, col=VISIT), jitter()) +
  geom_jitter(aes(colour=grouping))+
#  geom_boxplot(aes(fill=grouping))+
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +

  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw() 
  


#  stat_compare_means(comparisons = list(c("BL", "V08")), 
 #                    label = "p.format", method = "wilcox.test", tip.length = 0)






### compare and group
library(ggpubr)
library(rstatix)







ggplot(data = merged_melt_filt, aes(x = VISIT, y = value)) + 
  geom_point(aes(col=factor(VISIT)), size = 2) +
  geom_line(aes(group=patno, col=grouping)) +
  geom_boxplot(aes(fill=grouping))+
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  facet_wrap(. ~ symbol, scales='free_y') +
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw() + 
  stat_compare_means(comparisons = list(c("BL", "V08")), 
                     label = "p.format", method = "wilcox.test", paired=TRUE, tip.length = 0)






#dir.create(outdir, '/trajectories/')



#ggpaired(merged_melt, x = "VISIT", y = "value",
##        # color = "VISIT", line.color = "gray", line.size = 0.4,
 #        color = "VISIT", line.color = 'grouping', line.size = 0.4,
#         
#         palette = "jco")+
 # facet_wrap(facets = 'symbol', nrow=ceiling(length(ws_high)/3), ncol=4)+#
  #facet_free(facets = 'symbol')+
  
  #stat_compare_means(paired = TRUE)
  


merged_melt$patno


library(ggpubr)
# Create line plots of means
#ggline(merged_df, x = 'variable', y = "value", 
#       add = c("mean_sd", "jitter"))





