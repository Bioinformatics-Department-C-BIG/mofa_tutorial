




# Plot heatmap with dotted line for information of different groups 
# Factor X 
# Obtain for mofa plots 
# only factors we care about 
# 1. Variance captured: table
# 2. Cell type corelations log10 pvalue 
# 3. 

sel_factors_conf
estimations_in_df<-colnames(estimations)[colnames(estimations) %in% colnames(cors_all_pd)]
clinical_in_df<-c('NP2PTOT_LOG', 'NP3TOT_LOG',  'moca')
covars_age<-c('AGE_SCALED', 'SEX', 'LEDD')
sel_facts<-get_factors_for_scales(clinical_in_df)
#sel_facts<-sel_facts[!sel_facts %in% c(3)]
#sel_facts<-sel_facts[!sel_facts %in% c(3)]


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


library(circlize)
max_hm1<-max(abs(cors_cell_types))
max_hm2<-max(abs(cors_clinical))
max_hm3<-max(abs(cors_covars))
max_hm4<-max(abs(vars_factors))


col_fun1 = colorRamp2(c(0, max_hm1), c("white", "red"))
col_fun2 = colorRamp2(c(0, max_hm2), c("white", "orange"))
col_fun3 = colorRamp2(c(0, max_hm3), c("white", "orange"))

col_fun4 = colorRamp2(c(0, max_hm4), c("white", "purple"))

graphics.off()
jpeg(paste0(outdir, '/heatmap_factor_info.jpeg'), res=250, width=10, height=10, units='in')
cm1<-ComplexHeatmap::pheatmap(cors_cell_types, 
col =col_fun1)

cm2<-ComplexHeatmap::pheatmap(cors_clinical, col =col_fun2)
cm3<-ComplexHeatmap::pheatmap(cors_covars, col =col_fun3)

cm4<-ComplexHeatmap::pheatmap(vars_factors, col =col_fun4)


ht_list<- cm1 %v% cm2 %v% cm3 %v% cm4
draw(ht_list)

dev.off()








