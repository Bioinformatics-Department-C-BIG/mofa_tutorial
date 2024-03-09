

library(MOFA2)
library(dplyr)
library("vsn")


source(paste0('ppmi/setup_os.R'))
source(paste0('ppmi/config.R'));deseq_file;
source(paste0('ppmi/mofa_utils.R'))
source(paste0('ppmi/utils.R'))

output_files<-paste0(data_dir, '/ppmi/output/')

data_dir

# samples on columns 
metabolites<-read.csv(paste0(data_dir,'/thal/ThalCypris_peak_metaboanalyst_2022_corrected_all_sample_cells - Copy.csv'), row.names = 1)


mir_data<-read.csv(paste0(data_dir,'/thal/mature_counts.csv'), row.names = 1)

prot_data<-read.csv(paste0(data_dir,'thal/dia-proteinSummary-all_edited.csv'), row.names = 1)
colnames(prot_data)

samplesheet_original<-read.csv(paste0(data_dir,'thal/samples_metadata.csv'), row.names = 1)
#samplesheet<-process_clinvars(samplesheet)


### match the samplesheet to the first two? or the first two to the samplesheet? 




mir_data<-t(mir_data)


### TODO: RUN VSN ONLY make a standalone function 
colnames(metabolites)<-gsub(".*\\.","",colnames(metabolites))
colnames(prot_data); colnames(metabolites); colnames(mir_data)
samplesheet_original$Sample.name<-gsub(".*\\_","",samplesheet_original$Sample.name)



missing_sample<-samplesheet_original$Sample.name[!samplesheet_original$Sample.name %in% colnames(metabolites) ]

# Remove the missing sample - PS108
colnames(metabolites)
samplesheet<-samplesheet_original[!samplesheet_original$Sample.name %in% missing_sample,]



samplesheet$Sample.name

# 1. order metabolites by samplesheet , and match all others to metabolites - match them because they are the they are missing a sample 
metabolites_m<-metabolites[,na.omit(match(samplesheet$Sample.name, colnames(metabolites)))]
mir_data_m<-mir_data[,match(samplesheet$Sample.name, colnames(mir_data))]
prot_data_m<-prot_data[,match(samplesheet$Sample_prot, colnames(prot_data))]
colnames(prot_data_m)<-samplesheet$Sample.name
#colnames(mir_data_m)




colnames(mir_data)

samplesheet$sample<-samplesheet$sample.ID
samplesheet$COHORT<-as.factor(samplesheet$Biological.group)



#### 
dim(metabolites_m); colnames(metabolites_m)
dim(prot_data_m); colnames(prot_data_m)
dim(mir_data_m) ; colnames(mir_data_m)
metabolites2<-metabolites_m[-1,]

########## Preprocessing #### 

#metabolites2=metabolites_m
metabolites2<-metabolites2 %>% mutate_all(as.numeric)
prot_data_m<-prot_data_m %>% mutate_all(as.numeric)
### preprocessing
# TODO: functions 
metabolites_se<-SummarizedExperiment(metabolites2)

samplesheet$Biological.group
metabolites_se$COHORT<-samplesheet$Biological.group

mir_data_se<-SummarizedExperiment(mir_data_m)
mir_data_se$COHORT<-samplesheet$Biological.group


prot_data_se<-SummarizedExperiment(prot_data_m)
prot_data_se$COHORT<-samplesheet$Biological.group


########### MIRNAS preprocessing ######################

#### Do the filter by expression here! 
## TODO: check if this works now 
raw_counts=assays(mir_data_se)[[1]]
min.count=20
## filterbyExpr takes cpm so remove from there 
idx <- edgeR::filterByExpr(raw_counts,min.count=min.count, group = 'COHORT')
raw_counts <- as.matrix(raw_counts[idx, ])
dim(raw_counts)
mir_data_se=mir_data_se[idx]




formula_deseq = '~COHORT'
design = as.formula(formula_deseq)
ddsSE <- DESeqDataSet(mir_data_se,design )
ddsSE<-estimateSizeFactors(ddsSE)

vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
mirna_log<-log10(assay(ddsSE)+1)

### Metabolomics preprocessing 

metabolites_log<-as.matrix(log2(assay(metabolites_se)))



######### proteomics preprocessing #######
jpeg(paste0(data_dir, '/thal/plots/hist_proteins.jpeg'), units='in', height=10, width = 10, res=150)
hist(prot_data_m)
dev.off()
prot_data_m
prot_data_m<-as.matrix(prot_data_m)
colnames(prot_data_m)
#which(prot_data_m < 0)
prot_data_m <- replace(prot_data_m, which(prot_data_m < 0), NA)

jpeg(paste0(data_dir, '/thal/plots/hist_proteins_proc.jpeg'), units='in', height=10, width = 10, res=150)
hist(prot_data_m)
dev.off()


# dev.off()

#View(as.matrix(assay(vsd)))

### 
jpeg(paste0(data_dir, '/thal/plots/hist_proteins_metabolites.jpeg'), units='in', height=10, width = 10, res=150)
hist( as.matrix(metabolites_log))
dev.off()

jpeg(paste0(data_dir, '/thal/plots/hist_mirnas.jpeg'), units='in', height=10, width = 10, res=150)
hist( as.matrix(assay(vsd)))
dev.off()
jpeg(paste0(data_dir, '/thal/plots/hist_mirnas_log.jpeg'), units='in', height=10, width = 10, res=150)
hist( mirna_log)
dev.off()

####### HIGHLY VARIABLE######## 
most_var<-0.3
mirna_log_most_var =selectMostVariable(mirna_log, most_var)
metabolites_log_most_var =selectMostVariable(metabolites_log, most_var)
prot_data_m_most_var <- selectMostVariable(as.matrix(prot_data_m), most_var)

remove_outlier<-which(colnames(mirna_log_most_var)=='PS37')

######## run mofa ########################
data_full<-list(miRNA=as.matrix(mirna_log_most_var)[, -remove_outlier], 
                metabolites = as.matrix(metabolites_log_most_var)[, -remove_outlier], 
                proteomics=as.matrix(prot_data_m_most_var)[, -remove_outlier])


n=length(colnames(data_full[[1]]))

dim(mirna_log_most_var)

colnames(prot_data_m_most_var)
colnames(metabolites_log_most_var)

scale_views=TRUE

N_FACTORS=5
outdir<-paste0(data_dir,'/thal/output/three_mod_',n, '_',N_FACTORS, '/')

MOFAobject<-create_mofa(data_full)
MOFAobject<-run_mofa_wrapper( MOFAobject, outdir=outdir, N_FACTORS =N_FACTORS, force=FALSE)


MOFAobject
samplesheet_ro<-samplesheet[-remove_outlier,]
dim(samplesheet_ro)
samplesheet_ro$sample<-samplesheet_ro$Sample.name

samples_metadata(MOFAobject)<-samplesheet_ro















