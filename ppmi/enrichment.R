install.packages("msigdbr")
library(msigdbr)

subcategory='CP:KEGG';category='C2'
all_gene_sets = msigdbr(species = "human", category=category, subcategory = subcategory)


subcategory='GO:BP';category='C5'
all_gene_sets = msigdbr(species = "human", category=category, subcategory = subcategory)


colnames(all_gene_sets)
all_gene_sets_long<-all_gene_sets[,c('gs_name','ensembl_gene' )]

all_gene_sets_long<-all_gene_sets[,c('gs_name','ensembl_gene' )]
### better to obtain it from a protein db! 
all_gene_sets_long_prot<-all_gene_sets[,c('gs_name','gene_symbol' )]


print(all_gene_sets_long[, 'gs_name'], n=300)

all_unique_paths<-unique(all_gene_sets_long[, 'gs_name']) %>% as.data.frame()
colnames(all_unique_paths)
all_unique_paths$gs_name
all_unique_paths<-gsub('GOBP_', '', all_unique_paths$gs_name)
all_unique_paths
#Convert the VST counts to long format for ggplot2
library(reshape2)

library(pheatmap)
library(tidyr)
library(dplyr)

gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')

all_gene_sets_wide<-all_gene_sets_long %>%
  mutate(yesno = 1) %>%
  distinct %>%
  spread(ensembl_gene, yesno, fill = 0)

dim(all_gene_sets_wide)


all_gene_sets_wide


write.csv(all_gene_sets_wide,gs_file, row.names=FALSE )



#### proteins #### 
### created using gene symbols 

gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')
all_gene_sets_wide<-all_gene_sets_long_prot %>%
  mutate(yesno = 1) %>%
  distinct %>%
  spread(gene_symbol, yesno, fill = 0)
dim(all_gene_sets_wide)

write.csv(all_gene_sets_wide,gs_file, row.names=FALSE )

