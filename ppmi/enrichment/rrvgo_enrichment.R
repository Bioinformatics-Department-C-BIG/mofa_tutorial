
#BiocManager::install('rrvgo')
library(rrvgo)
ont='BP'
sem_sim<-GOSemSim::godata(orgdb, ont=ont)

go_analysis<-read.csv(paste0(data_dir,
  '/ppmi/plots/p_V08_CSF_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.2_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE/enrichment/gsego_12_.csv')
)


simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel", 
                                semdata = sem_sim
                                )


scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
# reduce the terms which are at least within a similarity belowe threshold
# and select the group representative 
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
