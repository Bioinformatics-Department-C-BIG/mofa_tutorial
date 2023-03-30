


output_1='ppmi/output/'
output_files_orig<-'ppmi/output/'
getwd()
disgenet_pd<-read.csv2('ppmi/ppmi_data/C0030567_disease_gda_summary_full.csv', sep = ',')

dir<-'ppmi/plots/p_BL-V04-V06-V08_Plasma_0.9_FALSE_1vsn_TRUE_g_BL-V04-V06-V08_0.1_100_1_m_BL-V04-V06-V08_0.5_10_1_8_coh_1_BL-V04-V06-V08/top_weights/'

dir<-'ppmi/plots/p_BL_Plasma_0.9_F_1vsn_TNA_0.8g_BL_0.1_100_1_m_BL_0.5_10_1_8_coh_1_BL_TRUE/top_weights/'

dir<-'ppmi/plots/p_BL_CSF_0.9_F_1vsn_TNA_0.8g_BL_0.1_100_1_m_BL_0.5_10_1_8_coh_1_BL_TRUE/top_weights/'

dir<-'ppmi/plots/p_BL_CSF_0.9_F_1vsn_TNA_0.8g_BL_0.1_100_1_m_BL_0.5_10_1_8_coh_1_BL_TRUE/top_weights/'

dir='ppmi/plots/p_BL_CSF_0.9_F_1-2vsn_TNA_0.8g_BL_0.1_100_1-2_m_BL_0.5_10_1-2_8_coh_1-2_BL_TRUE/top_weights/'

dir='ppmi/plots/p_BL_CSF_0.9_F_1-2vsn_TNA_0.8g_BL_0.1_100_1-2_m_BL_0.5_10_1-2_8_coh_1-2_BL_TRUE/top_weights/'

dir='2vsn_TNA_0.8g_BL_0.1_100_1'

files=list.files(dir, pattern = 'top', full.names = TRUE);files[5]
files2<-files[grep('miRNA', files, invert = TRUE)]

files2
T=0.01

truthset=disgenet_pd$Gene[disgenet_pd$Score_gda>0.01];truthset
file=files[5]
get_percent_indb<-function(file){
  
  
  #### 
      prot<-read.csv(file, sep='\t');prot
      top_feats<-prot$feature[prot$value>T]
      
      in_db<-top_feats %in% truthset
      percent_true<-mean(in_db)
      
      
      truthset_l=length(truthset)
      query_l=length(top_feats)
      
      overlap=length(which(in_db))
      popsize=length(unique(c(top_feats, truthset)))
     # print(paste(overlap, query_l, popsize- query_l,truthset_l))
      ph <-phyper(overlap,query_l,popsize- query_l,truthset_l, log.p=FALSE, lower.tail = FALSE )
      
      return(round(percent_true, digits=2))
      #return(ph)
      
}


res<-sapply(files2, get_percent_indb)
res
length(res)
mean(res, na.rm=TRUE)

#### OVEREPRESENTATION ANALYSIS 


res<-get_percent_indb(prot)
