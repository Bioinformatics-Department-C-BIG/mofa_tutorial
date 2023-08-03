

library(MOFA2)
library(dplyr)
source(paste0('ppmi/setup_os.R'))

source(paste0('ppmi/mofa_utils.R'))

source(paste0(script_dir, 'ppmi/config.R'));deseq_file;

# samples on columns 
metabolites<-read.csv('thal/ThalCypris_peak_metaboanalyst_2022_corrected_all_sample_cells - Copy.csv', row.names = 1)
metabolites

mir_data<-read.csv('thal/mature_counts.csv', row.names = 1)

prot_data<-read.csv('thal/dia-proteinSummary-all_edited.csv', row.names = 1)
colnames(prot_data)

samplesheet<-read.csv('thal/samples_metadata.csv', row.names = 1)
samplesheet<-process_clinvars(samplesheet)


### match the samplesheet to the first two? or the first two to the samplesheet? 




mir_data<-t(mir_data)

colnames(prot_data)
### TODO: RUN VSN ONLY make a standalone function 
colnames(metabolites)<-gsub(".*\\.","",colnames(metabolites))
samplesheet$Sample.name<-gsub(".*\\_","",samplesheet$Sample.name)

# 1. order metabolites by samplesheet , and match all others to metabolites - match them 
metabolites_m<-metabolites[,match(samplesheet$Sample.name, colnames(metabolites))]
mir_data_m<-mir_data[,match(samplesheet$Sample.name, colnames(mir_data))]
prot_data_m<-prot_data[,match(samplesheet$Sample_prot, colnames(prot_data))]
colnames(prot_data_m)<-samplesheet$Sample.name
#colnames(mir_data_m)




colnames(mir_data)

samplesheet$sample<-samplesheet$sample.ID
samplesheet$COHORT<-as.factor(samplesheet$Biological.group)



#### 

metabolites2<-metabolites_m[-1,]

########## Preprocessing #### 


metabolites2<-metabolites2 %>% mutate_all(as.numeric)
prot_data_m<-prot_data_m %>% mutate_all(as.numeric)
### preprocessing
# TODO: functions 
metabolites_se<-SummarizedExperiment(metabolites2)
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

### METAB preprocessing 

metabolites_log<-as.matrix(log2(assay(metabolites_se)))





########## proteomics preprocessing #######

hist(log10(as.matrix(prot_data_m)))


View(as.matrix(assay(vsd)))

### 
hist( as.matrix(metabolites_log))

hist( as.matrix(assay(vsd)))
hist( mirna_log)

####### HIGHLY VARIABLE######## 
most_var<-0.3
mirna_log_most_var =selectMostVariable(mirna_log, most_var)
metabolites_log_most_var =selectMostVariable(metabolites_log, most_var)

prot_data_m_most_var <- selectMostVariable(as.matrix(prot_data_m), most_var)



######## run mofa ########################
data_full<-list(miRNA=as.matrix(mirna_log_most_var), 
                metabolites = as.matrix(metabolites_log_most_var), 
                proteomics=as.matrix(prot_data_m_most_var))


dim(mirna_log_most_var)

colnames(prot_data_m_most_var)
colnames(metabolites_log_most_var)

scale_views=TRUE

outdir<-'./thal/output/three_mod/'
N_FACTORS=5
MOFAobject<-create_mofa(data_full)
MOFAobject<-run_mofa_wrapper( MOFAobject, outdir=outdir, N_FACTORS =N_FACTORS, force=TRUE)

samplesheet$sample<-samplesheet$Sample.name

samples_metadata(MOFAobject)<-samplesheet
MOFAobject
outdir













