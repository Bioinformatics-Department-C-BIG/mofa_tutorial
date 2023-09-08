


library(DescTools)




# Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model <- load_model(file)

# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 





clusters <- cluster_samples(MOFAobject, k=3, factors=c( 1,3,4,14))
clusters
#### TODO: tune with silhouette score ####



corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")



x1_all<-samples_metadata(MOFAobject)[selected_covars3][,1:10]

x2<-as.numeric(samples_metadata(MOFAobject)$cluster)
cor_clust <- psych::corr.test(samples_metadata(MOFAobject)[selected_covars2]$td_pigd_old_on ,as.numeric(samples_metadata(MOFAobject)$cluster), method = "pearson", adjust = "BH")
samples_metadata(MOFAobject)[selected_covars2]
#install.packages('DescTools')
apply(x1_all,2  ,MutInf, x2=x2, base=2)


for (i in 1: length(selected_covars3)){
  x1<-samples_metadata(MOFAobject)[selected_covars3][,i]
  print(paste(selected_covars3[i], MutInf(x1, x2)))
}




############## Create boxplots by group #### 


col_data<-samples_metadata(MOFAobject)[c(selected_covars2, 'cluster', 'PATNO')]
#col_data<-samples_metadata(MOFAobject)[c(selected_covars2)]

col_data_melt<-melt(col_data, id=c('PATNO', 'cluster'))

ggplot(col_data_melt)+ 
  geom_boxplot(aes(y=value, group=cluster))+
  facet_wrap(~variable,scales =  'free')
  

### Means by group 
library(dplyr)
group_by(col_data, cluster) %>%
  summarise(across(everything(), .f=list(mean=mean, sd=sd), na.rm=TRUE)
  
            
  )

group_by(col_data, cluster) %>%
  summarise(across(everything(), list(~var(., na.rm=TRUE)))
            
  )

 



