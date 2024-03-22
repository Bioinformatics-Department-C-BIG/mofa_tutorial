

source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_plots.R'))



# Plot heatmap with dotted line for information of different groups 
# Factor X 
# Obtain for mofa plots 
# only factors we care about 
# 1. Variance captured: table
# 2. Cell type corelations log10 pvalue 
# 3. 


estimations_in_df<-colnames(estimations)[colnames(estimations) %in% colnames(cors_all_pd)]
clinical_in_df<-c('NP2PTOT_LOG', 'NP3TOT_LOG',  'moca')
covars_age<-c('AGE_SCALED', 'SEX', 'LEDD')
sel_facts<-get_factors_for_scales(clinical_in_df)
#sel_facts<-sel_facts[!sel_facts %in% c(3)]
sel_facts<-sel_facts[!sel_facts %in% c(3)]


# Correlations: Clinical 
cors_cell_types<-t(cors_all_pd[sel_facts,estimations_in_df])

# Correlations: Cell types  
cors_clinical <-t(cors_all_pd[sel_facts,clinical_in_df])

# Correlations: Covariates 
cors_covars <-t(cors_all_pd[sel_facts,covars_age])




# Variances 
vars_factors<-t(vars_by_factor_all$r2_per_factor[[1]][sel_facts,])
vars_factors
length(cors_clinical)





########## Factor top weights ####

names(sel_facts)
select_top_bottom_perc(MOFAobject,view=1,   factors=21, top_fr=0.01)

# two options

MOFAobject
concatenate_top_features
top_ws_hm<-concatenate_top_features(MOFAobject, factors_all = names(sel_facts),view=1,top_fr = 0.01,weight_T = 0, pivot_wide = TRUE)
dim(top_ws_hm)
view = 'proteomics_csf'

top_ws_hm_prot<-concatenate_top_features(MOFAobject, factors_all = names(sel_facts),view=view,top_fr = 0.025,weight_T = 0.001, pivot_wide = TRUE)
names_vec =   rownames(top_ws_hm_prot)
dim(top_ws_hm_prot)


rownames(top_ws_hm_prot) =  modify_prot_names(names_vec, view, conv_uniprot = TRUE)
rownames(top_ws_hm_prot) 



library(circlize)
max_hm1<-max(abs(cors_cell_types))
max_hm2<-max(abs(cors_clinical))
max_hm3<-max(abs(cors_covars))
max_hm4<-max(abs(vars_factors))
max_hm5<-max(abs(top_ws_hm))


col_fun1 = colorRamp2(c(0, max_hm1), c("white", "red"))
col_fun2 = colorRamp2(c(0, max_hm2), c("white", "orange"))
col_fun3 = colorRamp2(c(0, max_hm3), c("white", "orange"))

col_fun4 = colorRamp2(c(0, max_hm4), c("white", "purple"))
col_fun5 = colorRamp2(c(0, max_hm5), c("white", "#16718f"))


cm1<-ComplexHeatmap::pheatmap(cors_cell_types, 
col =col_fun1, heatmap_legend_param = list(
        title='Correlation \n log10padj '
        ))
cors_cell_types

cm2<-ComplexHeatmap::pheatmap(cors_clinical, col =col_fun2, heatmap_legend_param = list(
        title='Correlation \n log10padj'
        ))
cm3<-ComplexHeatmap::pheatmap(cors_covars, col =col_fun3,heatmap_legend_param = list(
        title='Correlation \n log10padj'
        ))
cm4<-ComplexHeatmap::pheatmap(vars_factors, col =col_fun4, heatmap_legend_param = list(
        title='Var %'
))

top_ws_hm
cm5<-ComplexHeatmap::pheatmap(as.matrix(top_ws_hm),col =col_fun5, cluster_rows =  FALSE,
 heatmap_legend_param = list(
        title='Top features'
))

cm6<-ComplexHeatmap::pheatmap(as.matrix(top_ws_hm_prot),col =col_fun5, cluster_rows =  FALSE,
 heatmap_legend_param = list(
        title='Top features'
))



graphics.off()
jpeg(paste0(outdir, '/heatmap_factor_info.jpeg'), res=300, width=6, height=12, units='in')
ht_list<- cm1 %v% cm2 %v% cm3 %v% cm4%v% cm5%v% cm6
draw(ht_list)

dev.off()

















































































































































