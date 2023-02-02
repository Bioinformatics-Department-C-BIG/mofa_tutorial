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


outdir='ppmi/plots/'

### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')


# prerequisites: mass spec preprocessing and desq2 preprocessing

proteomics<-prot_bl_2

rownames(proteomics)<-proteomics$PATNO
proteomics$PATNO<-NULL


proteomics<-prot_bl_2_t
proteomics<-sapply(as.data.frame(proteomics),
                   as.numeric)

##### Handle duplicates 
proteomics<-proteomics[,!duplicated(colnames(proteomics),fromLast=TRUE)]

dim(proteomics)
length(unique(colnames(proteomics)))
vsd_mat<-selectMostVariable(proteomics,)



hist(as.numeric(proteomics))

miRNA<-highly_variable_genes_mofa



hist(as.numeric(miRNA))

dim(proteomics)

common_samples<-intersect(colnames(miRNA), colnames(proteomics))

which(colnames(miRNA) %in% common_samples)

miRNA_filt<-miRNA[,colnames(miRNA) %in% common_samples]
prot_filt<-proteomics[,colnames(proteomics) %in% common_samples]


dim(miRNA_filt)
dim(prot_filt)
## how to analyze and normalize olink? 



### INPUT TO MOFA
data = list(miRNA=log2(miRNA_filt), 
            proteomics = prot_filt)
#rownames(data$proteomics)<-paste(rownames(data$proteomics),"1",sep="_")



###  REMOVE NON ens ids 

order_cols<-colnames(miRNA_filt)
colnames(prot_filt)


select(prot_filt,order_cols)


hist(data$miRNA)
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

MOFAobject <- run_mofa(MOFAobject, outfile = 'mofa_bladder.hdf5')


pre_trained<-load_model('mofa_bladder.hdf5')
MOFAobject<-pre_trained
###

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


###
get_weights(MOFAobject)


