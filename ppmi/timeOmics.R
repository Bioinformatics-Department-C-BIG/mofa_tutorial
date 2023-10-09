
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
## install timeOmics
BiocManager::install('timeOmics')
BiocManager::install('propr')
#install.packages('propr')

library('timeOmics')
#data("timeOmics.simdata")

#install.packages("devtools")
# then load
library(devtools)
#install_github("abodein/timeOmics")
#remove.packages('timeOmics')
library('timeOmics')


mofa_multi_to_use

sim.data.se.all<-mofa_multi_to_use[ , ,'RNA']
sim.data.se.all<-mofa_multi_to_use[ , ,'miRNA']

order(colData(sim.data.se.all[c('PATNO', 'EVENT_ID'),] ))


# 




#combined_by_visit<-split(combined, combined$EVENT_ID )

met<-colData(sim.data.se.all)
byvisit<-split(met, met$EVENT_ID)
# add only 3 time points???? 
common_all_ts<-Reduce(intersect, x=list(byvisit$BL$PATNO,byvisit$V04$PATNO,byvisit$V06$PATNO,  byvisit$V08$PATNO ))
common_all_ts<-Reduce(intersect, x=list(byvisit$BL$PATNO,byvisit$V06$PATNO,  byvisit$V08$PATNO ))

length(common_all_ts)

sim.data.se.all.common<-sim.data.se.all[,colData(sim.data.se.all)$PATNO %in% common_all_ts]

sim.data.se.all.common<-sim.data.se.all.common[,colData(sim.data.se.all.common)$EVENT_ID %in% c('BL', 'V06', 'V08')]

sim.data.se.all.common<-sim.data.se.all.common[, order( colData(sim.data.se.all.common)[,c('PATNO', 'EVENT_ID')] )]
 head(colData(sim.data.se.all.common)[,c('PATNO', 'EVENT_ID')] )

 
sim.data.se.ct<-sim.data.se.all.common[,sim.data.se.all.common$COHORT==2]
sim.data.se.pd<-sim.data.se.all.common[,sim.data.se.all.common$COHORT==1]
 
sim.data=t(assay(sim.data.se))


# We assume each block (omics) is a matrix/data.frame with samples in rows (similar in each block) 
# and features in columns (variable number of column). 
remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X), 
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered.pd <- remove.low.cv(sim.data, 0.25)
data.filtered.ct <- remove.low.cv(t(assay(sim.data.se.ct)), 0.25)

colnames(data.filtered.ct)
colnames(data.filtered.pd)

common_temp_mols<-intersect(colnames(data.filtered.pd), colnames(data.filtered.ct))
data.filtered.only.pd <- data.filtered.pd[, !( colnames(data.filtered.pd) %in% common_temp_mols )  ]
data.filtered.only.ct<-data.filtered.ct[,!colnames(data.filtered.ct) %in% intersect(colnames(data.filtered.pd), colnames(data.filtered.ct))]

colnames(data.filtered.only.ct)



# remove the ones in pd
data.filtered=data.filtered.pd
sim.data.se=sim.data.se.pd

ens_ids<-colnames(data.filtered);
ens_ids
symbols_ids<-get_symbols_vector(ens_ids)
colnames(data.filtered)<-symbols_ids

dim(data.filtered )

data("timeOmics.simdata")
sim.data <- timeOmics.simdata$sim

dim(sim.data) 

sim.data

##### time modelling
#devtools::install_github("cran/lmms")
library(lmms)

# numeric vector containing the sample time point information


#### run by subtypes 

time <- factor(sim.data.se$EVENT_ID)

sim.data.se$PATNO


time<- factor(sim.data.se$EVENT_ID, levels=levels(as.factor(sim.data.se$EVENT_ID)),
              labels=c(seq(unique(sim.data.se$EVENT_ID))))

time<-as.numeric(time)
sid<-sim.data.se$PATNO_EVENT_ID
head(time)

any(is.na(data.filtered))
# example of lmms
lmms.output <- lmms::lmmSpline(data = data.filtered, time = time,
                               sampleID = sid, deri = FALSE,
                               basis = "p-spline", numCores = 1, timePredict = 1:max(time),
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))



# gather data
data.gathered <- modelled.data %>% as.data.frame() %>% 
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles




print(data.gathered[data.gathered$feature=='ADARB2', ])
p<-ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() +
  
  ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")
p
p + theme(legend.position = "none")


#Straight line modelling can occur when the inter-individual variation is too high. To remove the noisy profiles, 
#we have implemented a 2-phase test procedure.
# PROFILE FILTERING 
nrow(data.filtered)
length(time)
filter.res <- lmms.filter.lines(data = data.filtered, 
                                lmms.obj = lmms.output, time = as.numeric(time))
profile.filtered <- filter.res$filtered




# run pca
pca.res <- pca(X = profile.filtered, ncomp = 3, scale=FALSE, center=FALSE)

# tuning ncomp
pca.ncomp <- getNcomp(pca.res, max.ncomp = 3, X = profile.filtered, 
                      scale = FALSE, center=FALSE)

pca.ncomp$choice.ncomp

plot(pca.ncomp)




##### cluster
pca.res <- pca(X = profile.filtered, ncomp = 2, scale = FALSE, center=FALSE)
# extract cluster
pca.cluster <- getCluster(pca.res)
head(pca.cluster)

### multivar models 


plotIndiv(pca.res)
plotVar(pca.res)

plotLoadings(pca.res, comp = 1)


### PCA LONGITUDINAL 
plotLong(pca.res, scale = FALSE, center = FALSE, 
         title = "PCA longitudinal clustering")



############### LONGITUDINAL CLUSTERING 







