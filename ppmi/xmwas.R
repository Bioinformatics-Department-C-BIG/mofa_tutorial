



#### XMWAS ####
mirs_f<-paste0(output_files, '/xmwas/mirnas_V08_0.9_10_coh_1-2_INEXPD_highly_variable_genes_mofa.csv')
genes_f<-paste0(output_files, '/xmwas/rnas_V08_0.3_100_coh_1-2_INEXPD_highly_variable_genes_mofa.csv')
prot_f<-paste0(output_files, '/xmwas/V08_Plasma_0.3_T_1-2INEXPDvsn_TNA_0.9_highly_variable_proteins_mofa.csv')

outmirs<-paste0(output_files, '/xmwas/mirnas_V08_0.9_10_coh_1-2_INEXPD_highly_variable_genes_mofa_com.csv')
outgenes<-paste0(output_files, '/xmwas/rnas_V08_0.3_100_coh_1-2_INEXPD_highly_variable_genes_mofa_com.csv')
outgenes_prot<-paste0(output_files, '/xmwas/rnas_V08_0.3_100_coh_1-2_INEXPD_highly_variable_genes_mofa_com_prot.csv')


mirs<-read.csv(mirs_f, row.names = 1)
genes<-read.csv(genes_f, row.names = 1)
prot<-read.csv(prot_f, row.names = 1)


colnames(mirs)<-gsub('X', '', colnames(mirs))
colnames(genes)<-gsub('X', '', colnames(genes))
colnames(prot)<-gsub('X', '', colnames(prot))

com<-intersect(colnames(mirs), colnames(genes) )
com2<-intersect(colnames(prot), colnames(genes) )


length(mirs[, com])
length(genes[, com])
write.csv(mirs[, com],outmirs )
write.csv(genes[, com],outgenes )
write.csv(genes[, com2],outgenes_prot )

MET<-samples_metadata(MOFAobject)[,c('PATNO_EVENT_ID', 'COHORT')]
MET<-MET[MET$FileName%in%com,]
colnames(MET)<-c('FileName', 'Class')
write.csv(  MET, 
            file = paste0(output_files, '/xmwas/classes.csv'), 
            row.names = FALSE)




#####
ass_mat<-read.csv(paste0(output_files, '/xmwas/xmwasresults20230930143319/pairwise_results/datasetA_x_datasetB_association_matrix_threshold0.4.txt'), sep='\t')
ass_mat<-ass_mat[-1,-1]

ass_mat<-as.matrix(ass_mat)
library(igraph)
g <- graph.bipartite(ass_mat )


g <- graph.adjacency(ass_mat, )
get.edgelist(g)
