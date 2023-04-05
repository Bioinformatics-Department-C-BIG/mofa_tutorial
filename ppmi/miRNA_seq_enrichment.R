
#install.packages('rbioapi')
library(rbioapi)
order_by_metric<-'abslog2pval'
padj_T=0.01
log2fol_T=0.2

gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T, log2fol_T )
top_n=length(gene_list);
#top_n=300

gene_list_cut<-gene_list[1:top_n]
length(gene_list); length(gene_list_cut)
mirs=names(gene_list_cut)
mirs=names(gene_list)

mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                              mirna_type = "mature",
                              test_type = "GSEA",
                              species = 'Homo sapiens'
                              )

table(mieaa_all_gsea$Category)
mieaa_gsea_1<-mieaa_all_gsea[mieaa_all_gsea$Category=='Annotation (Gene Ontology)',]
mieaa_gsea_2<-mieaa_all_gsea[mieaa_all_gsea$Category=='GO Biological process (miRPathDB)',]
mieaa_gsea_3<-mieaa_all_gsea[mieaa_all_gsea$Category=='Gene Ontology (miRWalk)',]


dim(mieaa_gsea_1);
dim(mieaa_gsea_2);
dim(mieaa_gsea_3)

mieaa_gsea_2_cut<-mieaa_gsea_2[mieaa_gsea_2$`P-adjusted`<0.01, ]
dim(mieaa_gsea_2_cut)
mir_paths<-mieaa_gsea_2_cut[,c(2)]

results_file<-paste0(outdir_s, '/mirs_enrich_', '_', padj_T, '_',  log2fol_T, '_',  order_by_metric, '_',top_n)

write.csv(mieaa_all_gsea, paste0(results_file, '_', '.csv' ))



rna_p<-'ppmi/plots/single/rnas_V08_0.1_100_coh_1-2_AGE_AT_VISIT+SEX+COHORT/gseGO_BP_1_0log2pval.csv'
rna_paths<-read.csv(rna_p)
dim(rna_paths_cut)

rna_paths_cut<-rna_paths[rna_paths$p.adjust<0.01,]
common_pathways<-intersect(rna_paths_cut$Description ,mir_paths)
mir_paths[!(mir_paths %in% common_pathways)]
length(common_pathways)



