

# corelation of top gene markers with scores


# 

# todo: CHANGE THE VALUES of vsn to use RAW 



DIFF_VAR='NP3TOT_LOG'

DIFF_VAR='NP2PTOT_LOG'

DIFF_VAR='updrs_totscore'
DIFF_VAR='NP3TOT_LOG'
DIFF_VAR='sft'
DIFF_VAR='NP3TOT'




fact<-get_factors_for_metric(DIFF_VAR)

print(fact)

dir.create(paste0(outdir, '/correlations/'))



view='miRNA'
view='RNA'

view='proteomics_csf'


top_fr=0.1
MOFAobjectPD
top_genes<-select_top_bottom_perc(MOFAobjectPD,factor=fact[1], view=view)
length(top_genes)
top_genes_all<-concatenate_top_features(MOFAobjectPD,factors=fact, view=view, top_fr=top_fr)

top_genes<-top_genes_all$feature
top_genes_all[top_genes_all$feature == 'hsa-miR-101-3p', ]
head(get_data(MOFAobjectPD, view=view, groups=1, as.data.frame=TRUE))


patient_data<-as.data.frame(get_data(MOFAobjectPD, view=view, groups=1)[[1]])
patient_factor_data<-t(as.data.frame(MOFA2::get_factors(MOFAobjectPD)[[1]])[,fact])

dim(patient_factor_data)
patient_factor_data
colnames(patient_data)<-gsub('group1.','', colnames(patient_data) )
colnames(patient_data)
patient_factor_data
 patient_data = rbind(patient_data, patient_factor_data)




patient_data_top<-patient_data[rownames(patient_data) %in% c(top_genes),]

patient_data_top =  rbind(patient_data_top, patient_factor_data)
patient_scores<-MOFAobjectPD@samples_metadata[,DIFF_VAR]



patient_scores
length(patient_data_top)
x = patient_data_top[1,]

missing_data<-sapply(patient_data_top, function(x){all(is.na(x))})




patient_data_top_filt<-patient_data_top[,!missing_data]


patient_data_top_filt['Factor13',]
patient_scores_filt<-patient_scores[!missing_data]

ens_true<-startsWith(rownames(patient_data_top_filt)[1], 'ENS')

patient_data_top_filt


if (ens_true ){


        rownames(patient_data_top_filt)<-get_symbols_vector(rownames(patient_data_top_filt))
}else if (view == 'proteomics_csf'){
           rownames(patient_data_top_filt)<- get_symbol_from_uniprot(rownames(patient_data_top_filt))$SYMBOL



}
rownames(patient_data_top_filt)
patient_data_top_filt
cor_test<-corr.test(as.numeric(x),patient_scores)

as.numeric(x)
patient_scores
as.numeric(x)
Y=patient_scores
X=as.numeric(x)
fit = lm(Y~X+0)

summary(fit)
dim(patient_data_top_filt)
length(patient_scores_filt)

adjust ='fdr'
result_cor_1<-corr.test(t(patient_data_top_filt), patient_scores_filt,adjust = adjust)


result_cor<-data.frame(result_cor_1[c('r', 'p', 'p.adj')])
head(result_cor)




#patient_data_top_filt['ITGA2B',]

#corr.test(t(patient_data_top_filt['ITGA2B',]), patient_scores_filt)$r

pval_sel = 'p'

pCutoff = ifelse(pval_sel == 'p.adj', 0.1, 0.05)

max_pval<-max(-log10(result_cor[, pval_sel]))+0.5

result_cor$p

#ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 1),
xlim = c(-(max(abs(result_cor$r))+0.2), max(abs(result_cor$r))+0.2 )
max(abs(result_cor$r))



fact_labels<-grep('Factor',rownames(result_cor), value=TRUE)
fact_labels = NULL

p<-EnhancedVolcano(result_cor, 
        x='r', 
        y=pval_sel, 
        lab=rownames(result_cor), 
        selectLab = fact_labels,
        pCutoff = pCutoff, 
        FCcutoff = 0.1,
        xlim =xlim, 
        ylim= c(0, max_pval), 
        xlab = bquote('Pearson correlation'), 
        title = paste('Correlation', DIFF_VAR, '-', view  ),


        pointSize = 4.0,
        labSize = 7


)
show(p)
p
ggsave(paste0(outdir, '/correlations/cors_',DIFF_VAR,'_', view,top_fr,'pval',pval_sel, fact_labels[1], '.jpeg'), dpi=300)

















