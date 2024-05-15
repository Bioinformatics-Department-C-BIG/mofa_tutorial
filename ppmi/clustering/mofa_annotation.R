

#source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))
#source(paste0(script_dir, 'ppmi/mofa_analysis_plots.R'))
library('OmnipathR')

grn <-
    import_transcriptional_interactions(
        resources = c('DoRothEA')
    )
tfs<-unique(grn$source_genesymbol)
# Plot heatmap with dotted line for information of different groups 
# Factor X 
# Obtain for mofa plots 
# only factors we care about 
# 1. Variance captured: table
# 2. Cell type corelations log10 pvalue 
# 3. 


estimations_in_df<-colnames(estimations)[colnames(estimations) %in% colnames(cors_all_pd)]
clinical_in_df<-c('NP2PTOT_LOG', 'NP3TOT_LOG',  'moca', 'sft')
covars_age<-c('AGE_SCALED', 'SEX', 'LEDD', 'ab_asyn', 'tau_asyn', 'abeta', 'Neutrophil.Lymphocyte')
sel_facts<-get_factors_for_scales(clinical_in_df)
sel_facts
sel_facts<-get_factors_for_scales(clinical_in_df)
sel_facts2<-get_factors_for_metric('NP3TOT_LOG')

sel_facts<-sel_facts[sel_facts %in% sel_facts2]
#sel_facts
#sel_facts<-sel_facts[!sel_facts %in% c(3)]
#sel_facts<-sel_facts[!sel_facts %in% c(3)]


# Correlations: Clinical 
cors_cell_types<-t(cors_all_pd[sel_facts,estimations_in_df])
cors_cell_types_all<-t(cors[sel_facts,estimations_in_df])


cors_cell_types<-cors_cell_types[rowSums(cors_cell_types)>0,]


# Correlations: Cell types  
cors_clinical <-t(cors_all_pd[sel_facts,clinical_in_df])

# Correlations: Covariates 
cors_covars <-t(cors_all_pd[sel_facts,covars_age])




# Variances 
vars_factors<-t(vars_by_factor_all$r2_per_factor[[1]][sel_facts,])
vars_factors
length(cors_clinical)





########## Factor top weights ####


# two options


top_ws_hm<-concatenate_top_features(MOFAobject, factors_all = names(sel_facts),view=1,top_fr = 0.009,weight_T = 0.01, pivot_wide = TRUE)
dim(top_ws_hm)

# Only TFS 
only_tfs=FALSE

top_fr = ifelse(only_tfs, 0.02, 0.001 )
top_fr
top_ws_hm_rna<-concatenate_top_features(MOFAobject, factors_all = names(sel_facts),view=2,top_fr = top_fr,weight_T = 0.2, pivot_wide = TRUE)
dim(top_ws_hm_rna)
ens_orig<-rownames(top_ws_hm_rna)

symbs_orig<-get_symbols_vector(rownames(top_ws_hm_rna))

symbs<-symbs_orig
symbs[duplicated(symbs)]<-paste0(symbs[duplicated(symbs)], '_1')
symbs[duplicated(symbs)]<-paste0(symbs[duplicated(symbs)], '_1')
rownames(top_ws_hm_rna)<-symbs

dim(top_ws_hm_rna); colnames(top_ws_hm_rna)

if (only_tfs){
       top_ws_hm_rna<- top_ws_hm_rna[symbs %in% tfs,]

}
# which are TFs OR regulons? ? 





view = 'proteomics_csf'

top_ws_hm_prot<-concatenate_top_features(MOFAobject, factors_all = names(sel_facts),view=view,top_fr = 0.01,weight_T = 0.01, pivot_wide = TRUE)
rownames(top_ws_hm_prot) =  modify_prot_names(rownames(top_ws_hm_prot), view, conv_uniprot = TRUE)
dim(top_ws_hm_prot)

top_ws_hm_prot
sel_facts



#top_paths<-concatenate_top_pathways_factors(as.numeric(sel_facts),top_p=3, view=FALSE)
#top_paths$p.adjust = -log10(top_paths$p.adjust)
#top_paths_factors<-top_paths %>% pivot_wider(names_from = factor, values_from=p.adjust) %>% as.data.frame()
#rownames(top_paths_factors)<-top_paths_factors$Description;top_paths_factors$Description=NULL

#colnames(top_paths_factors)<-names(sel_facts); 

#top_paths_factors[is.na(top_paths_factors)]=0


### other pathways analysis 
sel_facts
view='RNA'
get_top_paths_matrix<-function(view){
        top_paths<-concatenate_top_pathways_factors(as.numeric(sel_facts),top_p=4, view=view, prefix=FALSE)
        top_paths$p.adjust = -log10(top_paths$p.adjust)
        top_paths
        top_paths=top_paths[!duplicated(top_paths$Description),]

        top_paths_factors<-top_paths %>% pivot_wider(names_from = factor, values_from=p.adjust) %>% as.data.frame()
        top_paths_factors=top_paths_factors[!duplicated(top_paths_factors$Description),]


        rownames(top_paths_factors)<-top_paths_factors$Description;top_paths_factors$Description=NULL
        colnames(top_paths_factors)

        sel_facts_m<-sel_facts[!sel_facts %in% colnames(top_paths_factors)]


        ### some factors are missing so add them to agree with the rest 
        empty_cols<-data.frame(matrix(rep(0,dim(top_paths_factors)[1] * length(sel_facts_m)), ncol=length(sel_facts_m) ))
        if (length(empty_cols)>0){
                empty_cols
                colnames(empty_cols)<-sel_facts_m
                top_paths_factors<-cbind(top_paths_factors, empty_cols)
        }


        top_paths_factors[is.na(top_paths_factors)]<-0
        colnames(top_paths_factors)<-names(sel_facts)[match(colnames(top_paths_factors), sel_facts)]
        top_paths_factors<-top_paths_factors[match(names(sel_facts),colnames(top_paths_factors))]
        return(top_paths_factors)



}

top_paths_factors_prot_csf<-get_top_paths_matrix(view='proteomics_csf')
top_paths_factors_prot_plasma<-get_top_paths_matrix(view='proteomics_plasma')
top_paths_factors_RNA<-get_top_paths_matrix(view='RNA')
top_paths_factors_prot_t_plasma<-get_top_paths_matrix(view='proteomics_t_plasma')

top_paths_factors_prot_t_csf<-get_top_paths_matrix(view='proteomics_t_csf')

top_paths_factors_RNA
top_paths_factors_RNA
#top_paths_factors_rna<-get_top_paths_matrix(view='RNA')


# TODO: protein enrichment analysis 

# Heatmap format settings ####

library(circlize)
max_hm1<-max(abs(cors_cell_types))
max_hm2<-max(abs(cors_clinical))
max_hm3<-max(abs(cors_covars))
max_hm4<-max(abs(vars_factors))
max_hm5<-max(abs(top_paths_factors_prot_csf))

vars_factors




col_fun1 = colorRamp2(c(0, max_hm1), c("white", "red"))
col_fun2 = colorRamp2(c(0, max_hm2), c("white", "orange"))
col_fun3 = colorRamp2(c(0, max_hm3), c("white", "orange"))

col_fun4 = colorRamp2(c(0, max_hm4), c("white", "purple"))
col_fun5 = colorRamp2(c(0, max_hm5), c("white", "#16718f"))
col_fun_paths = colorRamp2(c(0, max_hm5), c("white", "#168f48"))


cm1<-ComplexHeatmap::pheatmap(cors_cell_types, 
col =col_fun1, heatmap_legend_param = list(
        title='Correlation \n log10padj '
        ))
cors_cell_types

cm2<-ComplexHeatmap::pheatmap(cors_clinical, col =col_fun2, heatmap_legend_param = list(
        title='Correlation \n log10padj'
        ))
cm_covars<-ComplexHeatmap::pheatmap(cors_covars, col =col_fun3,heatmap_legend_param = list(
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

cm7<-ComplexHeatmap::pheatmap(as.matrix(top_ws_hm_rna),col =col_fun5, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='Top features'
))

cm_pathways<-ComplexHeatmap::pheatmap(as.matrix(top_paths_factors_prot_csf),col =col_fun_paths, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='csf proteomics'
))

cm_pathways_plasma<-ComplexHeatmap::pheatmap(as.matrix(top_paths_factors_prot_plasma),col =col_fun_paths, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='Plasma proteomics'
))

cm_rna<-ComplexHeatmap::pheatmap(as.matrix(top_paths_factors_RNA),col =col_fun_paths, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='RNA'
))

cm_t_csf<-ComplexHeatmap::pheatmap(as.matrix(top_paths_factors_prot_t_csf),col =col_fun_paths, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='targeted csf proteomics'
))


cm_t_plasma<-ComplexHeatmap::pheatmap(as.matrix(top_paths_factors_prot_t_plasma),col =col_fun_paths, cluster_rows =  FALSE, 
clustering_method = 'median',
 heatmap_legend_param = list(
        title='targeted plasma proteomics'
))
graphics.off()


jpeg(paste0(outdir, '/heatmap_factor_info',paste0(sel_facts, collapse='_'),only_tfs, '.jpeg'), res=300, width=6, height=6, units='in')
ht_list<-cm2  %v%  cm_covars %v% cm1  %v%  cm4   
draw(ht_list)

dev.off()


graphics.off()
jpeg(paste0(outdir, '/heatmap_factor_info_mol',only_tfs, '.jpeg'), res=300, width=10, height=17, units='in')
ht_list<-cm_rna %v% cm_pathways %v% cm_pathways_plasma %v% cm_t_csf %v% cm_t_plasma %v% cm7 %v% cm5%v% cm6 
draw(ht_list, padding = unit(c(2,2,2,80), 'mm'),  annotation_legend_side  = "bottom", 
heatmap_legend_side="bottom")


#draw(ch, heatmap_legend_side="bottom", 
 # annotation_legend_side  = "bottom", 
 #   padding = unit(c(2, 2, 2, 70), "mm"))
dev.off()



# TODO change the way we see pathways and add by factor? 





### More annotations 
#install.packages('hpar')
library('hpar')

library('hpar')
ens_orig1<-gsub('\\..*', '',ens_orig )
symbs<-get_symbols_vector(ens_orig)

import_omnipath_annotations(
  proteins = ens_orig1,
  resources = c('HPA_tissue'),
  wide = TRUE
  ) 

# import omnipath annotation : which tissue is the gene enriched? 
all_hpa<-import_omnipath_annotations(
  resources = c('HPA_tissue'),
  wide = TRUE
  ) 



colnames(all_hpa)
length(symbs)
gene_annotations<-all_hpa[all_hpa$genesymbol  %in% symbs,]
unique(gene_annotations$level)
gene_annotations<-gene_annotations[gene_annotations$status %in% c('Enhanced', 'Approved'),]

gene_annotations <-gene_annotations%>% 
                dplyr::filter(status %in% c('Enhanced', 'Approved') )  %>%
                dplyr::filter(!level %in% c('Not detected') ) 

gene_annotations <-gene_annotations%>% 
                dplyr::filter(status %in% c('Enhanced', 'Approved') )  %>%
                dplyr::filter(level %in% c('High', 'Medium') ) 

print(gene_annotations, n=100)
unique(gene_annotations$genesymbol)
table(gene_annotations$tissue)
colnames(gene_annotations)
print(gene_annotations,n=20)
gene_annotations
print(gene_annotations[grep('Purkinje|neur',gene_annotations$tissue ),], n=50)

gan_long<-gene_annotations[, c('genesymbol', 'tissue', 'level')]

gan_long2<-gan_long[!duplicated(gan_long),]
gan_long2$level[gan_long2$level=='High']=1
gan_long2$level[gan_long2$level=='Medium']=0.5
gan_long2$level<-as.numeric(gan_long2$level)

gan_long2<-gan_long2[!duplicated(gan_long2[,c('tissue', 'genesymbol')]),]




gan_long_wide<-gan_long2 %>% pivot_wider(names_from ='tissue', 
        values_from ='level') %>% as.data.frame()


rownames(gan_long_wide)<-gan_long_wide$genesymbol
gan_long_wide$genesymbol<-NULL
gan_long_wide[is.na(gan_long_wide)]=0
#apply(gan_long_wide, 2, as.numeric)

gan_long_wide
gan_long_wide[,grep('Purkinje',colnames(gan_long_wide))]
#gan_long_wide<-as.data.frame(gan_long_wide)
dim(gan_long_wide)


gan_long_wide2<-sapply(gan_long_wide,as.numeric)
rownames(gan_long_wide2)<-rownames(gan_long_wide)
rownames(gan_long_wide)


gan_long_wide2[,grep('Purkinje',colnames(gan_long_wide2))]


jpeg(paste0(outdir, '/gene_annotations.jpeg'), units='in', res=90, width=10,height=10)
pheatmap(t(gan_long_wide2), show_rownames = TRUE)
dev.off()


































































































































