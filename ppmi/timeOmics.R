
library('timeOmics')
data("timeOmics.simdata")

library('timeOmics')

#install.packages('cli')
library('cli')
sim.data <- timeOmics.simdata$sim

dim(sim.data) 



# load this from the full mofa run for VISIT=V08
#outdir_V08 = outdir
outdir_V08<-'/Volumes/GoogleDrive/Other computers/My computer (1) (1)/ppmi/plots/p_V08_CSF_T_1-2INEXPDvsn_TNA_0.9g_0.2_100_m_0.5_10_15_sig_0c_0_coh_1-2_V08_1ruv_1_c_1_split_0/'
all_clusts_V08<-read.csv2(paste0(outdir_V08,'/clustering/all_clusts_mofa.csv'), row.names=1, header=TRUE)

## run by cluster 
clust_ids<-all_clusts_V08[,'NP2PTOT_LOG']
names(clust_ids)<-rownames(all_clusts_V08)
clust_ids


# load this from dataset VISIT=all together  (partial run no full mofa run needed! )

mofa_multi_to_use_all_time<-mofa_multi_to_use




## 1. Run for all patients 
## 2. Run per group of patients 

mt<-colData(mofa_multi_to_use_all_time)
clust_ids_ordered<-clust_ids[match(mt$PATNO, names(clust_ids)  )]
mofa_multi_to_use_all_time$clust_ids<-clust_ids_ordered



clust_id_to_use=1
view='miRNA'
mofa_multi_to_use_all_time_clust<-mofa_multi_to_use_all_time[,mofa_multi_to_use_all_time$clust_ids %in% clust_id_to_use]
mofa_multi_to_use_all_time_clust


sim.data.se.all<-mofa_multi_to_use_all_time_clust[ , ,view]

#sim.data.se.all<-mofa_multi_to_use_all_time_clust[ , ,'proteomics_plasma']

assay(sim.data.se.all)
sim.data.se.all
sim.data.se<-sim.data.se.all[,sim.data.se.all$COHORT==1]
sim.data.se.all



sim.data=t(assay(sim.data.se))
head(sim.data)








#### WARNING:: DO NOT EDITTHE CODE BELOW. ediut THE SIM DATA FROM now on because everything is related to it ###
##

cutoff_rna<-0.1
cutoff_mirna<-0.22
# We assume each block (omics) is a matrix/data.frame with samples in rows (similar in each block) and features in columns (variable number of column). 
remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X), 
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered <- remove.low.cv(sim.data, cutoff_mirna)


dim(data.filtered)








##### time modelling
#devtools::install_github("cran/lmms")
library(lmms)

# numeric vector containing the sample time point information
time <- factor(sim.data.se$EVENT_ID)
time<- factor(sim.data.se$EVENT_ID, levels=levels(as.factor(sim.data.se$EVENT_ID)),
              labels=c("1", "2", "3", '4'))
time<-as.numeric(time)
sid<-sim.data.se$PATNO
head(time)
data.filtered2<-data.filtered

dim(data.filtered)
length(time)





lmms.output <- lmms::lmmSpline(data = data.filtered, time = time,
                               sampleID = sid, deri = FALSE,
                               basis = "p-spline", numCores = 1, timePredict = 1:4,
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))



# gather data
data.gathered <- modelled.data %>% as.data.frame() %>% 
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles
dev.off()
data.gathered$time
if(view=='RNA'){
  data.gathered$feature_ens<-get_symbols_vector(gsub('\\..*', '',data.gathered$feature))
}else{
  data.gathered$feature_ens = data.gathered$feature
}



p <- ggplot(data.gathered, aes(x = time, y = value, color = feature_ens)) + geom_line() +
  theme_bw() +
  
  ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")+
  theme(text=element_text(size=20))
p

#dir.create(paste0(outdir, '/time_omics/'))
ggsave(paste0(outdir, '/time_omics/',view, '_clust_', clust_id_to_use, '.png' ))


#Straight line modelling can occur when the inter-individual variation is too high. To remove the noisy profiles, 
#we have implemented a 2-phase test procedure.
# PROFILE FILTERING 
nrow(data.filtered)
length(time)
filter.res <- lmms.filter.lines(data = data.filtered, 
                                lmms.obj = lmms.output, time = as.numeric(time))
profile.filtered <- filter.res$filtered




# run pca
pca.res <- pca(X = profile.filtered, ncomp = 2, scale=FALSE, center=FALSE)

# tuning ncomp
pca.ncomp <- getNcomp(pca.res, max.ncomp = 2, X = profile.filtered, 
                      scale = FALSE, center=FALSE)

pca.ncomp$choice.ncomp

plot(pca.ncomp)

pca.ncomp


# TODO: continue time omics - time related for 
# 1. rna, mirs, proteomics
# run multi time omics 















