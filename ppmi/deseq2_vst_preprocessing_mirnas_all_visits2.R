

#### This script performs filter by expression, size factor estimation and vsn on the whole matrix of 
#### RNAs or miRNAS 
#### parameters: MIN_COUNT_M, MIN_COUNT_G define the parameters for filtering out genes with low counts 
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(sys)
library(DESeq2)
library("vsn")
library("SummarizedExperiment")
library(data.table)
library(dplyr)
## Output directory

output_1='ppmi/output/'
output_files_orig<-'ppmi/output/'
output_files<-'ppmi/output/'
outdir_orig<-'ppmi/plots/'

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)


output_de=paste0(output_1, 'gene')
source(paste0(script_dir, '/../bladder_cancer/preprocessing.R'))
source(paste0(script_dir, '/utils.R'))
source(paste0(script_dir, '/config.R'))


metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)

##### Load required data 
# TODO: input all the visits 

## TODO: move to config
## CONFIGURATION 


filter_common=TRUE

# MOVE ALL this to a configuration file!! 
#### Remove low expression 

process_mirnas<-TRUE
if (process_mirnas){
   mirnas_file<-paste0(output_files, 'mirnas_all_visits.csv')
   mirnas_BL<-as.matrix(fread(mirnas_file, header=TRUE), rownames=1)
   
   raw_counts<-mirnas_BL

  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, param_str_m, '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  deseq_file<-paste0(output_files, param_str_m,'deseq.Rds')
  
  
}else{
  rnas_file<-paste0(output_files, 'rnas_all_visits.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
 
  raw_counts<-rnas_BL
    # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
# raw_counts<-raw_counts %>% select(common_samples)

  min.count=MIN_COUNT_G
  most_var=TOP_GN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, 'rnas_', param_str_g,  '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  
  highly_variable_outfile
  
  
  deseq_file<-paste0(output_files, param_str_g,'deseq.Rds')
  
  
}
deseq_file
raw_counts_all=raw_counts


##### 1.  First create the summarized experiment object  
### find common samples in mirnas file + metadata
## subset and order by common samples
## And create SE object with metadata
#duplicate_samples<-colnames(mirnas_BL)[which(duplicated(colnames(mirnas_BL),fromLast = TRUE))]

class(raw_counts_all) <- "numeric"
## They seem to have taken averages for replicas so need to fix 
raw_counts_all<-round(raw_counts_all)

## Question: why are there duplicate samples - seems to be controls! 
## first filter what is in metadata and mirnas ?
se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined)
# remove duplicates 
##### Up till here it is generic, no filters yet. 


se_filt<-filter_se(se, VISIT, sel_coh)

### OUTPUT THE FILTERED se_filt 

ind<-which(is.na(se_filt$AGE_AT_VISIT))
se_filt[,ind]$AGE_AT_VISIT<-get_age_at_visit(colData(se_filt[,ind]))

## Turn to factors for deseq
se_filt$SEX<-as.factor(se_filt$SEX)


if (length(sel_coh)>1){
  
      if (length(VISIT)>1){
        print('Two cohorts and visits detected, running deseq and vsd with design formula')
        
        ddsSE <- DESeqDataSet(se_filt, 
                              design = ~COHORT + EVENT_ID  + SEX)
        ddsSE<-estimateSizeFactors(ddsSE)
        
        vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
        print(dim(vsd))
        
      }else{
        print('Two cohorts detected, running deseq and vsd with design formula')
      ddsSE <- DESeqDataSet(se_filt, 
                            design = ~SEX+COHORT )
      ddsSE<-estimateSizeFactors(ddsSE)
      
      vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
    
  }
  }else{
    print('Single cohort and visit deseq ')
    
    ddsSE <- DESeqDataSet(se_filt, 
                          design = ~PATNO + AGE_AT_VISIT + SEX)
    ddsSE<-estimateSizeFactors(ddsSE)
    
    vsd <- varianceStabilizingTransformation(ddsSE)
    
    
  }
deseq2Data <- DESeq(ddsSE)

datalist=list(ddsSE, vsd, se_filt,deseq2Data )
saveRDS(datalist,deseq_file)

se_filt$COHORT


# Compute normalization factors and vst 
# or use blind=false 

#ddsSE <- estimateSizeFactors(ddsSE)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function


vsd_mat <- assay(vsd)

###TODO: Move this to the mofa file 
highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
write.csv(highly_variable_genes_mofa, highly_variable_outfile, col.names = TRUE)
dim(highly_variable_genes_mofa)
rownames(highly_variable_genes_mofa)

highly_variable_outfile

run_plots<-FALSE
if (run_plots){
  meanSdPlot(vsd_mat)
  
  
  ######  Checks
  # TODO: could move checks outside pipeline in interactive mode
  # Check the effect of vst before and after
  par(mfrow=c(1,3))
  
  # Check distributions of samples using boxplots
  boxplot(log2(assay(ddsSE)), xlab="", ylab="Log2 counts ",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  title("Boxplots of logCPMs (unnormalised)")
  boxplot(log10(raw_counts), xlab="", ylab="Log10 counts ",las=2)
  
  # Check distributions of samples using boxplots
  boxplot(vsd_mat, xlab="", ylab="vst(counts) ",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(vsd_mat),col="blue")
  title("Boxplots of logCPMs (after vst)")
  
  ## ASSESS BATCH
  ### Assess batch effect
  #plotPCA(vsd, "sample") + labs(color='sample') + ggtitle("Batch effect") 
  #dev.off()
  
  #### This saves the whSole file without filtering for highly variable 
  write.csv(vsd_mat,vsn_out_file)
  ##### Store the most variable genes only for MOFA 
  # Select most variable genes
  ### run on its own for all visits? 
  # Check that the distribution is approximately normal
  dev.off()
  
  ### SANITY CHECK: Just plot one gene before and after preprocessing to ensure the mapping looks correct 
  df=highly_variable_genes_mofa
  par(mfrow=c(2,1))
  idx=30
  plot(df[idx,1:150]) 
  gname<-rownames(df)[idx]
  title(gname)
  df=raw_counts
  plot(df[gname, 1:150])
  title(gname)
  
}

