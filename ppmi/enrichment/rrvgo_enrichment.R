
#BiocManager::install('rrvgo')
library(rrvgo)
library('R.filesets')
ont='BP'
orgdb="org.Hs.eg.db"

##sem_sim<-GOSemSim::godata(orgdb, ont=ont)
#saveRDS(sem_sim,paste0('semsim', ont,'.RDS'))
# SAVE similarities to speed up reuse
sem_sim2<-loadRDS(paste0('semsim', ont,'.RDS'))
#sem_sim2<-loadRDS(paste0('semsim','.RDS'))

go_analysis<-read.csv(paste0(data_dir,
                             '/ppmi/plots/p_V08_CSF_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.2_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/gsego_12_.csv')
)

#pvalueCutoff_sig=0.05
cluster_id=3
results_file_cluster=paste0(cluster_params_dir, '/enr/'  ,'/gse',prefix, enrich_params, 'cl', cluster_id)


results_file_cluster_sig<-paste0(results_file_cluster, pvalueCutoff_sig, '.csv')
results_file_cluster_sig
go_analysis<-read.csv(results_file_cluster_sig)
go_analysis_sig<-go_analysis %>% dplyr::filter(p.adjust<pvalueCutoff_sig)
dim(go_analysis_sig)
print(go_analysis_sig)
# calculate sim matrix and reuse the sem_sim already calculate for this orddb and ont
simMatrix <- calculateSimMatrix(go_analysis_sig$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel", 
                                semdata = sem_sim2
                                )




scores <- setNames(-log10(go_analysis_sig$p.adjust), go_analysis_sig$ID)
length(scores)
# reduce the terms which are at least within a similarity belowe threshold
# and select the group representative 
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.6,
                                orgdb=orgdb)

dim(reducedTerms)
reducedTerms
graphics.off()
jpeg(paste0(results_file_cluster,'hm.png'))
p<-heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=1)


p
dev.off()
p<-scatterPlot( simMatrix,   reducedTerms, 
                     labelSize = 5       )+
        theme_gray()
p


ggsave(paste0(results_file_cluster,'scatterplot.png'))


