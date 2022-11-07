if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)


outdir='bladder_cancer/plots/mofa/'

### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')

mofa_proteins<-X2_raw

hist(highly_variable_proteins_mofa)
hist(highly_variable_genes_mofa)

data = list(mRNA=highly_variable_genes_mofa, 
            proteomics = highly_variable_proteins_mofa)
rownames(data$proteomics)<-paste(rownames(data$proteomics),"1",sep="_")



###  REMOVE NON ens ids 


head(colnames(reactomeGS))
symbols<-rownames(data$mRNA)

ens_ids <- mapIds(org.Hs.eg.db, keys =symbols, keytype = "SYMBOL", column="ENSEMBL")
rownames(data$mRNA)<-ens_ids
data$mRNA<-data$mRNA[-which(is.na(ens_ids)),]
###


hist(data$mRNA)
hist(data$proteomics)
data$proteomics
sapply(data$mRNA, as.numeric)

MOFAobject <- create_mofa(data)
MOFAobject



plot_data_overview(MOFAobject)



#### model opts 

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 4





model_opts
MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts
)




####Add metadata
bladder_metadata=Y_raw
bladder_metadata$sample<-samples_metadata(MOFAobject)$sample

samples_metadata(MOFAobject)<-bladder_metadata

MOFAobject <- run_mofa(MOFAobject)

saveRDS(MOFAobject, 'mofa.RDS')


###

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


###



