# Plot the heatmap
BiocManager::install('ComplexHeatmap')

library('ComplexHeatmap')

anno1<- c("purple","orange")[Y_raw$Subtype]
anno2 <- c("1","2")[as.factor(Y_raw$EORTC.risk)]
groups=Y_raw[c('EORTC.risk', 'Subtype', 'TURB.stage', 'Tumor.size')]

###### Load all genes and all features and clusters 

# best_clusters
#gene_cluster settings
ng_g=25
ng_p=4
k_genes=2;n.lambda_genes=89
k_proteins=2;n.lambda_proteins=13; upper_quantile=0.75
 #TODO: check if the names are correct ng, np

# Load clusters and features 
param_str_icl<-paste0(  most_var, '_ng_g_', ng_g,'_ng_p_', ng_p  )
fname_genes = paste0(icluster_out,'Vars_genes',param_str_icl, '_X_', k_genes, n.lambda_genes)
fname_proteins = paste0(icluster_out,'Vars_prot',param_str_icl, '_X_', k_proteins, n.lambda_proteins)
fname_proteins

k_multi=2; n.lambda_multi=233; ng_p=4; ng_g=25;upper_quantile=0.75
param_str_icl<-paste0(  most_var, '_ng_g_', ng_g,'_ng_p_', ng_p  )
fname_multi = paste0(icluster_out,'Vars_multi',param_str_icl, '_X_', k_multi, n.lambda_multi)

gene_clusters= data.frame(read.csv(paste0(fname_genes,  'best_clusters.csv' ),row.names=1))[,1]
protein_clusters=read.csv( paste0(fname_proteins,  'best_clusters.csv' ), row.names = 1)[,1]
multi_clusters=read.csv( paste0( fname_multi,'best_clusters.csv'  ) , row.names = 1)[,1]
multi_clusters_CIMLR=  c(1, 1, 2, 1, 2, 2, 3, 1, 3, 3, 2, 3, 3, 3, 3, 3)

top_proteins=read.csv(paste0(fname_proteins, '_', upper_quantile,  '.csv' ),row.names=1)[,1]
top_genes= read.csv(paste0(fname_genes, '_', upper_quantile,  '.csv' ),row.names=1)[,1]



column_ha = groups






ha=HeatmapAnnotation(EORTC.risk=Y_raw$EORTC.risk,
                   Subtype=Y_raw$Subtype,
                   TURB.stage= Y_raw$TURB.stage,
                   Grade=Y_raw$Grade, 
                   #Size=as.factor(Y_raw$Tumor.size),
                   iClusterGene=c(gene_clusters),
                   iClusterProt=c(protein_clusters),
                   iClusterMulti=c(multi_clusters),
                   cimlrMulti=c(multi_clusters_CIMLR),
                   
                     gp = gpar(col = "black"))
highly_variable_genes_voom_m<-data.matrix(highly_variable_genes_voom)
highly_variable_proteins_voom_m<-data.matrix(highly_variable_proteins_voom)

cluster_rows=TRUE
png(paste0(icluster_out, 'heatmap_all_', cluster_rows, '_',  k_multi, '_',upper_quantile, '.png'), height=900, width=500, res = 100)

Heatmap(highly_variable_genes_voom_m[top_genes,], 
        cluster_columns = TRUE,
        cluster_rows = cluster_rows,
          top_annotation = ha, 
        row_names_gp = gpar(fontsize = 7))

dev.off()
top_proteins %in% rownames(highly_variable_proteins_voom_m)
Heatmap(highly_variable_proteins_voom_m[top_proteins,], top_annotation = ha)




### new heatmap using both features from cilmr 
png(paste0(icluster_out, 'heatmap_all_cimlr_', cluster_rows, '_',  k_multi, '_',upper_quantile, '.png'), height=900, width=500, res = 100)

in_genes<-which(rownames(input_data)[head(ranks$aggR,50)] %in% rownames(highly_variable_genes_voom_m))
in_proteins<-which(rownames(input_data)[head(ranks$aggR,50)] %in% rownames(highly_variable_proteins_voom_m))
rownames(highly_variable_proteins_voom_m)[rownames(input_data)[head(ranks$aggR,40)] %in% rownames(highly_variable_proteins_voom_m)]

Heatmap(rbind(highly_variable_genes_voom_m[in_genes,], highly_variable_proteins_voom_m[in_proteins,]), 
        cluster_columns = TRUE,
        cluster_rows = cluster_rows,
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 7))

dev.off()

