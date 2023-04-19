

#### This script performs filter by expression, size factor estimation and vsn on the whole matrix of 
#### RNAs or miRNAS 
#### parameters: MIN_COUNT_M, MIN_COUNT_G define the parameters for filtering out genes with low counts 

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(sys)
library(GenomicRanges)

#remove.packages('Glimma') 
library(DESeq2)
library("SummarizedExperiment")
library(data.table)
library(dplyr)
## Output directory
output_de=paste0(output_1, 'gene')


script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))

source(paste0(script_dir, '/../bladder_cancer/preprocessing.R'))
source(paste0(script_dir, '/utils.R'))

metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)




#for (VISIT in c('V08', 'BL')){
  for (VISIT in c('V08')){
    
  
  source(paste0(script_dir, '/config.R'))
  
  raw_counts<-as.matrix(fread(input_file, header=TRUE), rownames=1)
  raw_counts_all<-raw_counts
  
  ##### Load required data 
  # TODO: input all the visits 
  
  ## TODO: move to config
  ## CONFIGURATION 
  
  
  filter_common=TRUE
  
  # MOVE ALL this to a configuration file!! 
  #### Remove low expression 
  
  
  
  
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
  
  
  ### I moved this elsewhere 
  # TODO: use se_filt$AGE_SCALED and test!!
  se_filt<-filter_se(se, VISIT, sel_coh)
  ### OUTPUT THE FILTERED se_filt 
  
  ind<-which(is.na(se_filt$AGE_AT_VISIT))
  se_filt[,ind]$AGE_AT_VISIT<-get_age_at_visit(colData(se_filt[,ind]))
  
  ## Turn to factors for deseq
  se_filt$SEX<-as.factor(se_filt$SEX)
  
  se_filt$AGE_AT_VISIT<-scale(se_filt$AGE_AT_VISIT)
  
  
  
  if (length(sel_coh)>1){
    
    if (length(VISIT)>1){
      print('Two cohorts and visits detected, running deseq and vsd with design formula')
      
      ddsSE <- DESeqDataSet(se_filt, 
                            design = as.formula(formula_deseq))
      ddsSE<-estimateSizeFactors(ddsSE)
      
      vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
      
      
      
    }else{
      print('Two cohorts detected, running deseq and vsd with design formula')
      ddsSE <- DESeqDataSet(se_filt, 
                            design =as.formula(formula_deseq2 ))
      ddsSE<-estimateSizeFactors(ddsSE)
      
      
      ### separate vsd? 
     # se_filt[]
     # vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
      vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
      
      
    }
  }else{
    print('Single cohort and visit deseq ')
    
    ddsSE <- DESeqDataSet(se_filt, 
                          design = formula_deseq3)
    ddsSE<-estimateSizeFactors(ddsSE)
    
    vsd <- varianceStabilizingTransformation(ddsSE)
    
    
  }
  
  deseq2Data <- DESeq(ddsSE)
  deseq2Results <- results(deseq2Data, contrast=c('COHORT', 1,2))
  datalist=list(ddsSE, vsd, se_filt, deseq2Results)
  saveRDS(datalist,deseq_file)
  
  
  
  
  # Compute normalization factors and vst 
  ### select mofa genes 
  # or use blind=false 
  
  # ddsSE <- estimateSizeFactors(ddsSE)
  # Variance stabilization transformation
  # This uses the size factors estimated before 
  # TODO: you can run VST using a saved dispersion function
  
  # in mofa apply also a filtering based on most DE genes
  #
  
  
  ## need to move to deseq analaysi 
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
  
  
  
  
  
}

