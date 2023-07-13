
library(org.Hs.eg.db)
library(edgeR)
source(paste0(script_dir, 'ppmi/utils.R'))


## Load baseline and v08 ####
# run deseqvst
bl_f<-'ppmi/plots/p_BL_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_BL_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
bl<-read.csv(bl_f)
bl<-bl[!is.na(bl$Description),]
bl_sig<-bl[bl$Least_value<0.01,]

v08_f<-'ppmi/plots/p_V08_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
v08<-read.csv(v08_f)


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
view='proteomics'; process_mirnas=FALSE ## NEED TO LOAD proteins df for this -- TODO: fix olink preprocesing
view='miRNA'; process_mirnas=TRUE
view='RNA'; process_mirnas=FALSE

source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
se_filt_V06<-filter_se(se, VISIT='V06', sel_coh,sel_ps)
se_filt_V04<-filter_se(se, VISIT='V04', sel_coh,sel_ps)

#Reduce(intersect, list(a,b,c))

common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )




ws<-get_weights(MOFAobject, views = view, factors=sel_factors[fn_sel])[[1]]


if (fn_sel==2){
  cut_high<-0.98; cut_low=0.02
  
}else{
  cut_high<-0.998; cut_low=0.002
  
}
ws_high<-ws[ws>quantile(ws, cut_high),]
ws_low<-ws[ws<quantile(ws, cut_low),]
ws_high

ws_union<-c(ws_high, ws_low)
#ws_high<-ws_high[1:15]
length(ws_union)
#ws_high<-ws_high[abs(ws_high)<quantile(ws_high, 0.9)]






if (view=='RNA'){
  ens_genes<-rownames(assay(se_filt_V08))[grep(paste0(names(ws_union), collapse='|'), rownames(assay(se_filt_V08)))]
  feat_names=ens_genes
  
}else{
  feat_names=names(ws_union)
  
}

function(){
  
}





preprocess_visit<-function(se_filt_V){
  
  se_filt_V_pd<-se_filt_V[,se_filt_V$COHORT == 1]
  se_filt_V_pd<-se_filt_V_pd[,se_filt_V_pd$PATNO %in% common]
  df_v<-cpm(assay(se_filt_V_pd),  normalized.lib.sizes=TRUE, log=TRUE )
  df_v<- clip_outliers(df_v)
  df_V_ens<-t(df_v[rownames(df_v) %in% feat_names,])
  v_ens=data.frame(df_V_ens,patno=c(se_filt_V_pd$PATNO), PATNO_EVENT_ID=c(se_filt_V_pd$PATNO_EVENT_ID))
  
  return(v_ens)
  
}

v6_ens<-preprocess_visit(se_filt_V06)
v8_ens<-preprocess_visit(se_filt_V08)
se_filt_V08_pd<-se_filt_V08[,se_filt_V08$COHORT == 1]


v4_ens<-preprocess_visit(se_filt_V04)
bl_ens<-preprocess_visit(se_filt_BL)







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








######### PLOT molecular markers 



### MELT and MERGE 
v8_melt<-melt(v8_ens)
v6_melt<-melt(v6_ens)
v4_melt<-melt(v4_ens)
bl_melt<-melt(bl_ens)


bl_melt$VISIT<-'BL'
v4_melt$VISIT<-'V04'
v8_melt$VISIT<-'V08'
v6_melt$VISIT<-'V06'



merged_melt<-rbind(bl_melt, v4_melt)
merged_melt<-rbind(merged_melt,v6_melt)
merged_melt<-rbind(merged_melt,v8_melt)

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
Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,0.8, na.rm=TRUE))

sel_factors[fn_sel]
if (names(sel_factors[fn_sel]) %in% c('Factor4', 'Factor14', 'Factor1')){
  T=0.2
  Z1_grouping<-factor(Z1_matched>quantile(Z1_matched,T, na.rm=TRUE))
  
}
T
Z1_grouping
dim(merged_melt)
merged_melt$grouping<-Z1_grouping
merged_melt

if (view=='RNA'){
  symb<-get_symbols_vector(ens)
  merged_melt$symbol<-symb
}else{
  merged_melt$symbol<-merged_melt$variable
  
}

na_ps<-unique(merged_melt[!is.na(merged_melt$grouping),]$patno)
merged_melt_filt<-merged_melt[merged_melt$patno %in% na_ps, ]

### EDIT GROUP TO BE PATNO and GROUP!!
# OR JUST GROUP  
merged_melt_filt$VISIT


group_up<-unique(names(Z1_grouping[which(Z1_grouping==FALSE)]))
group_down<-unique(names(Z1_grouping[which(Z1_grouping==TRUE)]))
### Did they change scales???

se_filt_V08_pd_g1<-se_filt_V08_pd[,se_filt_V08_pd$PATNO_EVENT_ID %in% group_down ]
se_filt_V08_pd_g2<-se_filt_V08_pd[,se_filt_V08_pd$PATNO_EVENT_ID %in% group_up ]


G1<-se_filt_V08_pd_g1[,!(se_filt_V08_pd_g1$PDSTATE == 'ON')]
G2<-se_filt_V08_pd_g2[,!(se_filt_V08_pd_g2$PDSTATE == 'ON')]
## TODO: plot here ALL the clinical variables by grouping!! 

#### 

merged_melt_filt$VISIT<-as.factor(merged_melt_filt$VISIT)

merged_melt_filt$grouping<-as.factor(merged_melt_filt$grouping)









################


### Plot to remove the other group ####
# TAKE THE low group
if (names(sel_factors[fn_sel]) %in% c('Factor4', 'Factor14', 'Factor1')){
  GROUP=FALSE
  
  
}else{
  GROUP=TRUE
}
merged_melt_filt_g1=merged_melt_filt[merged_melt_filt$grouping %in% c(GROUP),]
merged_melt_filt_g1=merged_melt_filt_g1[merged_melt_filt_g1$VISIT %in% c('BL', 'V08'),]
merged_melt_filt_g2=merged_melt_filt[merged_melt_filt$grouping %in% c(GROUP),]

merged_melt_filt_g1$VISIT<-as.factor(merged_melt_filt_g1$VISIT)
merged_melt_filt_g2$VISIT<-as.factor(merged_melt_filt_g2$VISIT)



wilcox_stats<-merged_melt_filt_g1 %>% group_by(symbol) %>%
  do(w=wilcox.test(value~VISIT, data=.))%>%
  summarize(symbol, Wilcox=w$p.value) %>%
  as.data.frame()

most_sig_over_time<-wilcox_stats[order(wilcox_stats$Wilcox),][1:5,]

merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time$symbol,]

ggplot(data = merged_melt_filt_g2_sig, aes(x = VISIT, y = value)) + 
  geom_point(aes(col=VISIT), size = 2) +
  geom_line(aes(group=patno),  col= 'grey') +
  geom_boxplot(aes(fill=VISIT))+
  scale_color_viridis_d(option='mako')+
  scale_fill_viridis_d(option='mako')+
  
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  
  geom_signif(comparisons = list(c('BL', 'V08')),  
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.4)+

  facet_wrap(. ~ symbol, scales='free_y', nrow=1) +
  
    theme_bw() 
ggsave(paste0(outdir, '/trajectories/boxplots_', sel_factors[fn_sel],'_', view,'_', GROUP, '.jpeg'), 
       width=12, height=3)



#######################################################
############ TIME TRAJECTORY FOR ALL VISITS #####



mean_data<-merged_melt_filt_g1 %>%
  group_by(grouping, symbol) %>%
  summarise(mean_expr = mean(value, na.rm = TRUE))
mean_data
  most_sig_over_time
  
  median_IQR <- function(x) {
    data.frame(y = median(x), # Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
  }
  

  filt_top=TRUE
if (filt_top){
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% most_sig_over_time$symbol[1:5],]
}else{
  merged_melt_filt_most_sig<-merged_melt_filt
  
}

ggplot(data = merged_melt_filt_most_sig, aes(x = VISIT, y = value, 
                                    fill=grouping, group=grouping, colour=grouping)) + 
  stat_summary(geom = "pointrange", fun.data = median_IQR, 
               position=position_dodge(0))+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1) + 
  scale_color_viridis_d(option='turbo')+
  facet_wrap(. ~ symbol, scales='free_y', 
            nrow = 1) +
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw()+ 
 geom_signif(comparisons = list(c('BL', 'V08')), 
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.4)+
  
  labs(y='logCPM')+
 # legend(legend=c('Low', 'High'))+
  theme(strip.text = element_text(
    size = 12, color = "dark green"), 
    axis.title.y =element_text(
      size = 12, color = "dark green") )
  
  

warnings()
ggsave(paste0(outdir, '/trajectories/trajectory_', sel_factors[fn_sel],'_', view, filt_top,   '.jpeg'), 
       width=10, height=3)












stat.test2 <- merged_melt_filt %>%
  group_by(grouping) %>%
  t_test(value ~ VISIT, p.adjust.method = "bonferroni")
stat.test2

##############

### Plot clinical variables ####

na_ps
### AND ALSO changed STAGE 

se_V08_cl<-se_filt_V08_pd[,se_filt_V08_pd$PATNO%in% na_ps ]
se_bl_cl<-se_filt_BL_pd[,se_filt_BL_pd$PATNO%in% na_ps ]


se_V08_cl$NP4TOT
sel_factors[fn_sel]
cors_both<-get_correlations(MOFAobject, names(MOFAobject@samples_metadata))
cors_all=cors_both[[1]]

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
#table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'],merged_melt_cl$grouping )
#table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'], )
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
merged_melt_cl3<-merged_melt_cl[merged_melt_cl$PDSTATE %in% c('OFF', ''),]
merged_melt_cl2<-merged_melt_cl3
merged_melt_cl3$PDSTATE

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





