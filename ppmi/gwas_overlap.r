
data_dir
gwas_genes<-read.csv(paste0(script_dir,'/ppmi/replication/gwas_genes.csv'), header=FALSE)


fact_cohort<-which(cors[, 'COHORT']>0);fact_cohort

diff_var='sft'
fact_cohort<-get_factors_for_metric('sft')
fact_cohort<-get_factors_for_metric('NP2PTOT_LOG')
fact_cohort<-get_factors_for_metric('sft')
fact_cohort<-which(cors[, 'COHORT']>0);fact_cohort

fact_cohort<-get_factors_for_metric(diff_var)
fact_cohort<-which(cors[, 'COHORT']>0);fact_cohort
fact_cohort
view = 'proteomics_csf'; top_fr=0.1
view = 'RNA' ; top_fr = 0.05

view = 'proteomics_t_plasma'; top_fr=0.1

top_rna_feats<-concatenate_top_features(MOFAobject,factors_all=fact_cohort, top_fr = top_fr, view=view)

top_rna_feats$feature<-convert_to_gene_symbol(top_rna_feats$feature, view=view)
top_rna_feats$feature  # factor 3
overlapping<-intersect(gwas_genes$V1,top_rna_feats$feature)
print(overlapping)
top_rna_feats[top_rna_feats$feature %in% overlapping,]



