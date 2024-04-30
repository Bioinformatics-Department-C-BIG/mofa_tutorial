

# corelation of top gene markers with scores


# 


DIFF_VAR='NP3TOT_LOG'

DIFF_VAR='NP2PTOT_LOG'

fact<-get_factors_for_metric(DIFF_VAR)

fact
view='miRNA'
view='proteomics_csf'
view='RNA'

top_fr=0.05
MOFAobjectPD
top_genes<-select_top_bottom_perc(MOFAobjectPD,factor=fact[1], view=view)
length(top_genes)
top_genes_all<-concatenate_top_features(MOFAobjectPD,factors=fact, view=view, top_fr=top_fr)

top_genes<-top_genes_all$feature


patient_data<-as.data.frame(get_data(MOFAobjectPD, view=view)[[1]])

patient_data
patient_data_top<-patient_data[rownames(patient_data) %in% top_genes,]
dim(patient_data_top)
patient_scores<-MOFAobjectPD@samples_metadata[,DIFF_VAR]



patient_scores
length(patient_data_top)
x = patient_data_top[1,]
patient_data_top
missing_data<-sapply(patient_data_top, function(x){all(is.na(x))})




patient_data_top_filt<-patient_data_top[,!missing_data]

patient_scores_filt<-patient_scores[!missing_data]

ens_true<-startsWith(rownames(patient_data_top_filt)[1], 'ENS')

if (ens_true ){


rownames(patient_data_top_filt)<-get_symbols_vector(rownames(patient_data_top_filt))
}

patient_data_top_filt
cor_test<-corr.test(as.numeric(x),patient_scores)
dim(patient_data_top_filt)
length(patient_scores_filt)
result_cor_1<-corr.test(t(patient_data_top_filt), patient_scores_filt,adjust = 'fdr')
head(result_cor)

result_cor<-data.frame(result_cor_1[c('r', 'p', 'p.adj')])



patient_data_top_filt['ITGA2B',]
patient_scores_filt
corr.test(t(patient_data_top_filt['ITGA2B',]), patient_scores_filt)$r

pval_sel = 'p'

max_pval<-max(-log10(result_cor[, pval_sel]))+0.1

result_cor$p

#ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 1),
xlim = c(-(max(abs(result_cor$r))+0.2), max(abs(result_cor$r))+0.2 )
max(abs(result_cor$r))

p<-EnhancedVolcano(result_cor, 
x='r', 
y=pval_sel, 
lab=rownames(result_cor), 
pCutoff = 0.05, 
FCcutoff = 0.2,
xlim =xlim, 
ylim= c(0, max_pval), 
xlab = bquote(correlation)

)

show(p)

ggsave(paste0(outdir, '/correlations_',DIFF_VAR,'_', view,top_fr, '.jpeg'), dpi=300)

















