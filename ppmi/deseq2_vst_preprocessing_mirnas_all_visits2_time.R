

#### This script performs filter by expression, size factor estimation and vsn on the whole matrix of 
#### RNAs or miRNAS 
#### parameters: MIN_COUNT_M, MIN_COUNT_G define the parameters for filtering out genes with low counts 
### The parameters used are defined in the script config.R

library(edgeR)
library(limma)
#library(Glimma)
#library(gplots)
library(RColorBrewer)
library(sys)
library(GenomicRanges)
#BiocManager::install('DeSeq2')
library(DESeq2)
library("SummarizedExperiment")
library(data.table)
library(dplyr)
library(rbioapi)


## Output directory
# output_de=paste0(output_1, 'gene')
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, '/bladder_cancer/preprocessing.R'))
source(paste0(script_dir, 'ppmi/utils.R'))

### Load metadata 
metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)
metadata_output<-paste0(output_files, 'combined_log.csv')
combined_bl_log<-read.csv2(metadata_output)

process_mirnas=TRUE
### Perform deseq for each visit (timepoint separately)
#for (VISIT in c('V08', 'BL')){

# tTODOl: FOR BL MATCH THE V08 SAMPLES!! DO
VISIT='V08'
VISIT=c('BL','V08')
VISIT=c('BL','V04', 'V06',  'V08');


        print(VISIT)
        filter_common=TRUE
        same_samples_all_visits=TRUE
        source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
        filter_by_group=TRUE
        
        ##### Load required data 
        # TODO: input all the visits 
        # MOVE ALL this to a configuration file!! 
        #### Remove low expression 
        
        
        se=load_se_all_visits(input_file = input_file, combined=combined_bl_log)
        
        ##### 1.  First create the summarized experiment object  
        ### find common samples in mirnas file + metadata
        ## subset and order by common samples
        ## And create SE object with metadata
        # remove duplicates 
        ##### Up till here it is generic, no filters yet. 
        se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
        se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
        
        common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )
       
        ### I moved the age scaling elsewhere 
        # TODO: use se_filt$AGE_SCALED and test!!
        se_filt<-filter_se(se, VISIT, sel_coh)
        se_filt<-filter_se(se, VISIT, sel_coh, sel_ps)
        
        se_filt_V08<-filter_se(se, VISIT='V08', sel_coh, sel_ps)
        se_filt_V04<-filter_se(se, VISIT='V04', sel_coh, sel_ps)
        se_filt_V06<-filter_se(se, VISIT='V06', sel_coh, sel_ps)
        se_filt_BL<-filter_se(se, VISIT='BL', sel_coh, sel_ps)
        se_filt_V04<-se_filt_V04[,!(is.na(se_filt_V04$SEX))]
        
        
        dim(se_filt)

        common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )
        if (same_samples_all_visits){
          common=Reduce(intersect, list(se_filt_V08$PATNO,se_filt_BL$PATNO , se_filt_V04$PATNO, se_filt_V06$PATNO))
          
        }
        
        
        
        
        if (filter_common){
          se_filt<-se_filt[,se_filt$PATNO %in% common]
        }
        dim(se_filt)
        #### Do the filter by expression here! 
        ## TODO: check if this works now 
        
       
        
        se_filt<-filter_se_byExpr(se_filt)
        
        ### TODO: ADD FILTER COMMON AS PARAM TO SAVE IN THE DIRECTORY IN CONFIG!!! 
        
        
        
        ### OUTPUT THE FILTERED se_filt 
        ind<-which(is.na(se_filt$AGE_AT_VISIT))
        se_filt[,ind]$AGE_AT_VISIT<-get_age_at_visit(colData(se_filt[,ind]))
        
        ## Turn to factors for deseq
        se_filt$SEX<-as.factor(se_filt$SEX)
        se_filt$AGE_AT_VISIT<-scale(se_filt$AGE_AT_VISIT)
        
        ## these are almost the same so it is okay to scale AGE earlier 
        hist(se_filt$AGE_AT_VISIT)
        hist(se_filt$AGE_SCALED)
        
        # impute: 
        # which()
        se_filt$AGE_SCALED[is.na(se_filt$AGE_SCALED)]<-mean(se_filt$AGE_SCALED, na.rm=TRUE)
        se_filt<-se_filt[,!(is.na(se_filt$SEX))]
        
        table(colData(se_filt)[,c( 'EVENT_ID', 'SEX')])

        colData(se_filt)[,c( 'EVENT_ID', 'SEX', 'AGE', 'PATNO')]
        
        ### Perform the appropriate test depending on what you want as prediction variable
        if (length(sel_coh)>1){
          
          if (length(VISIT)>1){
            print('Two cohorts and visits detected, running deseq and vsd with design formula')
             
            #se_filt2 = se_filt[, se_filt$COHORT==2]
            se_filt2 = se_filt
            
            if (same_samples_all_visits){
              se_filt2 = se_filt[, se_filt$COHORT==2]
              se_filt2=se_filt
             # gsub('\\_.*', '', names(kmeans_grouping))
              filter_by_group=FALSE
              if (filter_by_group){
                patnos_g1<-gsub('\\_.*', '', names(which(groups_kmeans$cluster==2)))
                
                cns<-se_filt[,se_filt$COHORT_DEFINITION=='Healthy Control']$PATNO
                se_filt2 = se_filt[, se_filt$PATNO %in% c(patnos_g1, cns)]
                
                
                
              }
              
              
              formula_deseq = '~AGE_SCALED+SEX+COHORT + EVENT_ID + COHORT:EVENT_ID'

             
              
            }
            
            ddsSE <- DESeqDataSet(se_filt2, 
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
                                design = as.formula('~AGE_AT_VISIT+SEX'))
          ddsSE<-estimateSizeFactors(ddsSE)
          
          vsd <- varianceStabilizingTransformation(ddsSE)
          
          
        }
        
        
      
        
        deseq2Data <- DESeq(ddsSE)
        ### Contrast disease-control: parkinsons = 1, control = 2 
        se_filt$COHORT_DEFINITION; se_filt$COHORT
        if (4 %in% sel_coh){
           disease_group=4    
        }else{
          disease_group=1    
          
        }
        
        if (length(sel_coh)>1){
          # Check if there is more than one cohort
          deseq2Results <- results(deseq2Data, contrast=c('COHORT', disease_group,2))
        }else{
          deseq2Results <- results(deseq2Data)
          
        }
        
        datalist=list(ddsSE, vsd, se_filt, deseq2Results)
        saveRDS(datalist,deseq_file)
        saveRDS(se, se_file)
        
        deseq2Results_sex <- results(deseq2Data, contrast=c('SEX', 0,1))
        
        deseq2Results_sex[order(deseq2Results_sex$padj),]
        
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
  
  
  
  
  


