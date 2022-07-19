# Plot the heatmap
#BiocManager::install('ComplexHeatmap')

library('ComplexHeatmap')

anno1<- c("purple","orange")[Y_raw$Subtype]
anno2 <- c("1","2")[as.factor(Y_raw$EORTC.risk)]
groups=Y_raw[c('EORTC.risk', 'Subtype', 'TURB.stage')]



# best_clusters
#gene_cluster settings
ng_g=35
ng_p=5
k_genes=2;n.lambda_genes=89
k_proteins=2;n.lambda_proteins=13
param_str_icl<-paste0(  most_var, '_ng_g_', ng_g,'_ng_p_', ng_p  )
fname_genes = paste0(icluster_out,'Vars_genes',param_str_icl, '_X_', k_genes, n.lambda_genes, 'best_clusters.csv')
fname_proteins = paste0(icluster_out,'Vars_prot',param_str_icl, '_X_', k_genes, n.lambda_proteins, 'best_clusters.csv')
fname_proteins

k_multi=3; n.lambda_multi=233; ng_p=4; ng_g=25
param_str_icl<-paste0(  most_var, '_ng_g_', ng_g,'_ng_p_', ng_p  )

fname_multi = paste0(icluster_out,'Vars_multi',param_str_icl, '_X_', k_multi, n.lambda_multi, 'best_clusters.csv')

gene_clusters= data.frame(read.csv(fname_genes,row.names=1))[,1]
protein_clusters=read.csv(fname_proteins, row.names = 1)[,1]
multi_clusters=read.csv(fname_multi, row.names = 1)[,1]

top_n=80
column_ha = groups

ha=HeatmapAnnotation(EORTC.risk=Y_raw$EORTC.risk,
                   Subtype=Y_raw$Subtype,
                   TURB.stage= Y_raw$TURB.stage,
                   Grade=Y_raw$Grade,
                   iClusterGene=c(gene_clusters),
                   iClusterProt=c(protein_clusters),
                   iClusterMulti=c(multi_clusters),
                   
                     gp = gpar(col = "black"))
highly_variable_genes_voom_m<-data.matrix(highly_variable_genes_voom)
highly_variable_proteins_voom_m<-data.matrix(highly_variable_proteins_voom)

Heatmap(highly_variable_genes_voom_m[top_genes[1:top_n],], 
        cluster_columns = TRUE,
          top_annotation = ha, 
        row_names_gp = gpar(fontsize = 6))

top_proteins %in% rownames(highly_variable_proteins_voom_m)
Heatmap(highly_variable_proteins_voom_m[top_proteins,], top_annotation = ha)


