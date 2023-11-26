
#BiocManager::install('rrvgo')
library(rrvgo)
library('R.filesets')
ont='BP'
orgdb="org.Hs.eg.db"

#sem_sim<-GOSemSim::godata(orgdb, ont=ont)
#saveRDS(sem_sim,paste0('semsim', ont,'.RDS'))
# SAVE similarities to speed up reuse
sem_sim2<-loadRDS(paste0('semsim', ont,'.RDS'))
sem_sim2<-loadRDS(paste0('semsim','.RDS'))

go_analysis<-read.csv(paste0(data_dir,
                             '/ppmi/plots/p_V08_CSF_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.2_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/gsego_12_.csv')
)

# calculate sim matrix and reuse the sem_sim already calculate for this orddb and ont
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel", 
                                semdata = sem_sim2
                                )

go_analysis$p.adjust
scores <- setNames(-log10(go_analysis$p.adjust), go_analysis$ID)
length(scores)
# reduce the terms which are at least within a similarity belowe threshold
# and select the group representative 
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=1.3,
                                orgdb="org.Hs.eg.db")

dim(reducedTerms)
p<-heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=1)


p

p<-scatterPlot(simMatrix, reducedTerms)

ggsave(paste0(outdir,'/enrichment/scatterplot.png'))
