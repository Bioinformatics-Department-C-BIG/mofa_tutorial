

#### another method:: 
#BiocManager::install('multiMiR')
library('multiMiR')
### get all possible targets from mirtar base

#BiocManager::install("targetscan.Hs.eg.db")
library(targetscan.Hs.eg.db)
f=4
for (f in 1:10){
  anti_rnas<-get_weights(MOFAobject, views='RNA', factors=f)[[1]]
  anti_rnas_2=anti_rnas[order(anti_rnas),][1:500]
  anti_rnas_2
  mirs_to_test<-get_weights(MOFAobject, views='miRNA', factors = f)[[1]]
  mirs_to_test_2=mirs_to_test[order(-mirs_to_test),][1:100]
  mirs_to_test_2
  
  ### check if mofa targets have the anticorelated rnas inside 
  library('org.Hs.eg.db')
  all_tars<-mget(names(mirs_to_test), revmap(targetscan.Hs.egTARGETS))
  ids_ens<-mapIds(org.Hs.eg.db, keys = unlist(all_tars), keytype="ENTREZID", column = "ENSEMBL")
  all_ens<-unlist(ids_ens)
  #example1 <- get_multimir(mirna =names(mirs_to_test_2), summary = TRUE)
  
  checks<-all_ens %in% c(anti_rnas_2)
  print(which(checks))
}




sexample1 <- get_multimir(mirna =mirs[N], summary = TRUE)
write.csv(example1@data, paste0(mir_results_file, 'gene_targets.csv'))


table(example1@data$type)
head(example1@data)

all_targets<-example1@data[c('mature_mirna_id', 'target_ensembl')]
anticor_long<-melt(anticor)

dim(all_targets)
colnames(all_targets)

anticor_long_ints<-anticor_long[anticor_long$value,]

colnames(anticor_long_ints)<-c('target_ensembl', 'mature_mirna_id', 'int')
anticor_long_ints$mature_mirna_id<-gsub('\\.', '-', anticor_long_ints$mature_mirna_id)

head(all_targets);head(anticor_long_ints)

### now filter all targets by what is actually anticorrelated 
merged_targets<-merge(all_targets, anticor_long_ints, by=c('mature_mirna_id', 'target_ensembl'))

write.csv(merged_targets, paste0(mir_results_file, 'gene_targets_filtered', N, '.csv'))

