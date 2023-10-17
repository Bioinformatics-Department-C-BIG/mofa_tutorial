
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
## install timeOmics
BiocManager::install('timeOmics')
BiocManager::install('propr')
install.packages('propr')

library('timeOmics')
data("timeOmics.simdata")

install.packages("devtools")
# then load
library(devtools)
install_github("abodein/timeOmics")
library('timeOmics')

install.packages('cli')
library('cli')
sim.data <- timeOmics.simdata$sim

dim(sim.data) 

mofa_multi_to_use


## 1. Run for all patients 
## 2. Run per group of patients 

sim.data.se.all<-mofa_multi_to_use[ , ,'miRNA']
sim.data.se.all<-mofa_multi_to_use[ , ,'RNA']

sim.data.se<-sim.data.se.all[,sim.data.se.all$COHORT==1]

sim.data=t(assay(sim.data.se))
head(sim.data)



# We assume each block (omics) is a matrix/data.frame with samples in rows (similar in each block) and features in columns (variable number of column). 
remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X), 
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered <- remove.low.cv(sim.data, 0.15)




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
data.filtered2<-na.omit(data.filtered)

any(is.na(data.filtered))
# example of lmms
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
p<-ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() +
  
  ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")
p
#p + theme(legend.position = "none")


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

