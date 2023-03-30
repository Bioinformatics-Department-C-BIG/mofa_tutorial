#### Comparisons:
# TODO compare: which are common 
# TODO: compare p values
# Now also compare V08 and CONTROLS ! 
# SOS: check if the comparisons are between the same patients !! 

#install.packages('venneuler')
#install.packages('sets')

output_comp<-'ppmi/plots/comparisons/'

library('venneuler')
library(sets)
V06_neg<-"ppmi/plots/p_V06_Plasma_0.9_g_V08_0.1_100_m_V06_0.5_10_10/GO_BP_enrichment_negative_pvals_no_f_p_V06_Plasma_0.9_g_V06_0.1_100_m_V06_0.5_10_10.csv"
V06_pos<-"ppmi/plots/p_V06_Plasma_0.9_g_V06_0.1_100_m_V06_0.5_10_10/GO_BP_enrichment_positive_pvals_no_f_p_V06_Plasma_0.9_g_V06_0.1_100_m_V06_0.5_10_10.csv"

neg<-read.csv(V06_neg)
pos<-read.csv(V06_pos)

V06_paths<-rbind(neg, pos)

V08_neg<-'ppmi/plots/p_V08_Plasma_0.9_g_V08_0.1_100_m_V08_0.5_10_10/GO_BP_enrichment_negative_pvals_no_f_p_V08_Plasma_0.9_g_V08_0.1_100_m_V08_0.5_10_10.csv'
V08_pos<-'ppmi/plots/p_V08_Plasma_0.9_g_V08_0.1_100_m_V08_0.5_10_10/GO_BP_enrichment_positive_pvals_no_f_p_V08_Plasma_0.9_g_V08_0.1_100_m_V08_0.5_10_10.csv'


# make sure we have the same patients!



neg<-read.csv(V08_neg)
pos<-read.csv(V08_pos)

V08_paths<-rbind(neg, pos)



BL_neg<-"ppmi/plots/p_BL_Plasma_0.9_g_BL_0.1_100_m_BL_0.5_10_10/GO_BP_enrichment_negative_pvals_no_f_p_BL_Plasma_0.9_g_BL_0.1_100_m_BL_0.5_10_10.csv"
BL_pos<-"ppmi/plots/p_BL_Plasma_0.9_g_BL_0.1_100_m_BL_0.5_10_10/GO_BP_enrichment_positive_pvals_no_f_p_BL_Plasma_0.9_g_BL_0.1_100_m_BL_0.5_10_10.csv"
neg<-read.csv(BL_neg)
pos<-read.csv(BL_pos)

BL_paths<-rbind(neg, pos)

colnames(V08_paths)[3]<-'paths'
colnames(BL_paths)[3]<-'paths'

NROW(unique(V08_paths$paths))
NROW(unique(BL_paths$paths))


common_paths_V08_BL<-intersect(V08_paths$paths,BL_paths$paths) # WHICH ARE Common without checking pvalues, just common paths 
common_paths_V06_BL<-intersect(V06_paths$paths,BL_paths$paths) # WHICH ARE Common without checking pvalues, just common paths 

NROW(unique(file1))
NROW(file1)

V08_only<-V08_paths[!(V08_paths$paths %in%  common_paths_V08_BL),]
V06_only<-unique(V06_paths[!(V06_paths$paths %in% common_paths_V06_BL),])
BL_only<-unique(file2[!(file2$paths %in% common_paths),])

write.csv(BL_only, paste0(output_comp,'BL_only.csv'))
write.csv(V08_only, paste0(output_comp,'V08_only.csv'))
write.csv(V06_only, paste0(output_comp,'V06_only.csv'))



