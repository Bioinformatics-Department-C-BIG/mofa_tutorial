

# corelation of top gene markers with scores


# 


DIFF_VAR='NP3TOT_LOG'

DIFF_VAR='NP2PTOT_LOG'

fact<-get_factors_for_metric(DIFF_VAR)

fact
view='miRNA'
view='proteomics_csf'
view='RNA'

top_fr=0.2
MOFAobjectPD
top_genes<-select_top_bottom_perc(MOFAobjectPD,factor=fact[1], view=view)
top_genes_all<-concatenate_top_features(MOFAobjectPD,factors=fact, view=view, top_fr=top_fr)

top_genes<-top_genes_all$feature
top_genes
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

cor_test<-cor.test(as.numeric(x),patient_scores)
cor_test$estimate
result_cor<-apply(patient_data_top_filt,1, function(x){

    cor_test<-cor.test(x,patient_scores_filt)
    return(c(cor_test$estimate, cor_test$p.value))
    })
result_cor<-as.data.frame(t(result_cor))
result_cor


max_pval<-max(-log10(result_cor$V2))+1



#ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 1),
xlim = c(-(max(abs(result_cor$cor))+0.2), max(abs(result_cor$cor))+0.2 )
max(abs(result_cor$cor))

p<-EnhancedVolcano(result_cor, 
x='cor', 
y='V2', 
lab=rownames(result_cor), 
pCutoff = 0.05, 
FCcutoff = 0.1,
xlim =xlim, 
ylim= c(0, max_pval), 
xlab = bquote(correlation)

)

show(p)

ggsave(paste0(outdir, '/correlations_',DIFF_VAR,'_', view,top_fr, '.jpeg'), dpi=300)















