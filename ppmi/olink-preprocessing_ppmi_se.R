
#script_dir='/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/'
#os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'

source('ppmi/setup_os.R')

data_dir

#data_dir<-os_dir
#BiocManager::install('DEP')
## TODO: change all scripts to be agnostic of visit until mofa
library(limma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(ggplot2)
library("vsn")
library("data.table")
library("SummarizedExperiment")

#### TODO: SAVE SE FILT SO WE CAN RELOAD in the next script 

#if (!require("pacman")) install.packages("pacman")
source(paste0(script_dir,'ppmi/utils.R'))
output_1=paste0(data_dir,'ppmi/plots/proteomics/')
outdir_orig<-paste0(data_dir,'ppmi/plots/')
output_files<-paste0(data_dir,'ppmi/output/')


metadata_output<-paste0(output_files, 'combined_log.csv')
combined<-read.csv2(metadata_output)
combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)





metadata_output<-paste0(output_files, 'combined_log.csv')
combined<-read.csv2(metadata_output)


### TODO: filter out the cohort too before processig !! 

VISIT=c('BL')
#process_mirnas=FALSE
#TISSUE='CSF'

source(paste0(script_dir, 'ppmi/config.R' ))

param_str<-paste0(TOP_PN)

if (NORMALIZED){
  in_file_original<-paste0(output_files, 'proteomics_', p_params_in,  '_no_log.csv')
  # if we dont run vsn we want to take the log values 
  if (!run_vsn){
    in_file_original<-paste0(output_files, 'proteomics_', p_params_in,  '.csv')
  }
}else{
  in_file_original<-paste0(output_files, 'proteomics_', p_params_in, '.csv')
}
in_file_original


#### Read in 
prot_bl_wide_unlog<-as.matrix(fread(in_file_original, header=TRUE), rownames=1)
in_file_original
proteomics<-prot_bl_wide_unlog
  
#proteomics[2,]
## TODO: tidy up!!

colnames(proteomics)


#### FILTERING LOW VALUES 
# Remove rows with 90% NA 

#### TODO: move this after the selection of the cohort!!! 

pre_process_proteomics<-function(proteomics){

      df<-proteomics
      proteomics <- df[rowSums(is.na(df)) < round(0.2*ncol(df)), ]

      
            ### filter here before editing more 
      ## filter out rows with very low min count
      df<-proteomics; 
      min.count= quantile(df, na.rm = TRUE, 0.01)
      min.count= min(df, na.rm = TRUE)
      ### KEEP
       ## kEEP THE ROWS THAT HAVE MORE THAN 80% NON NA VALUES 
      keep <- rowSums(df>min.count, na.rm = TRUE) >= round(NA_PERCENT*ncol(df))
      proteomics<-proteomics[keep,]

      
      
      #### MAKE NUMERIC 
      raw_counts_all=proteomics
      class(raw_counts_all) <- "numeric"
      ## They seem to have taken averages for replicas so need to fix 
      #raw_counts_all<-round(raw_counts_all)
      
      data<-proteomics
      data<-as.data.frame(data)
      
      data$name<-c(rownames(data))
      data$ID<-data$name
      data_columns=seq(1:dim(proteomics)[2])
      
      return(raw_counts_all)
}
#sample<-colnames(proteomics_se)


raw_counts_all<-pre_process_proteomics(proteomics)
proteomics_se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined)
dim(proteomics_se)

##### filter here by visits and cohort

se_filt<- filter_se(proteomics_se, VISIT, sel_coh, sel_ps)
se_filt_proteins<- filter_se(proteomics_se, VISIT, sel_coh, sel_ps)

### TODO: save se filt here : with or without VISIT included..? 
Sample<-colnames(se_filt)
sample_info<-DataFrame(Sample=Sample)

assays
tmp<- assays(se_filt)[[1]]


# Filter and normalize
### ERROR: Error in vsnML(sv) : L-BFGS-B needs finite values of 'fn'

### TODO: filter before normalization!!! 

normalized_data<-justvsn(tmp)
vsn::meanSdPlot(normalized_data)



ggsave(paste0(output_1,'meansd_justvsn_', p_params_out,'.png' ), width = 5, height=3)


vsn_mat<-normalized_data
head(vsn_mat)
# Select the top most variable proteins
## TODO: fix the bug in selectMostVariable

p_params_out<- paste0(VISIT_S, '_',TISSUE, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)


write.csv2(vsn_mat,paste0(output_files,p_params_out, '_vsn.csv'), row.names=TRUE)


for (most_var in c(0.3, 0.9)){
    highly_variable_proteins_outfile<-paste0(output_files, p_params_out , '_highly_variable_proteins_mofa.csv')
    highly_variable_sign_proteins_outfile<-paste0(output_files, p_params_out , '_highly_variable_proteins_mofa_signif.csv')

    highly_variable_proteins_mofa=selectMostVariable(vsn_mat, most_var)
}
png(paste0(output_1,'hist_high_var_', p_params_out,'.png' ))
hist(highly_variable_proteins_mofa)
dev.off()


png(paste0(output_1,'hist_', p_params_out,'.png' ))
hist(vsn_mat)
dev.off()



#### save all and load 
##deseq2Results <- results(deseq2Data, contrast=c('COHORT', 1,2))
datalist=list( vsn_mat, se_filt) # save the vsn summarized experiment and the raw summarized experiment 
prot_vsn_se_filt_file
saveRDS(datalist,prot_vsn_se_filt_file)
meanSdPlot(vsn_mat)



# IF WE DONT RUN VSN WE SAVE THE SE VALUES 
if (!run_vsn){
  print('Saving without VSN')
  highly_variable_proteins_mofa=selectMostVariable(tmp, TOP_PN)
  write.csv(highly_variable_proteins_mofa,highly_variable_proteins_outfile)
  meanSdPlot(tmp)
  ggsave(paste0(output_1,'meansd_NO_VSN_', p_params_out,'.png' ), width = 5, height=3)
  png(paste0(output_1,'hist_', p_params_out,'.png' ))
  hist(tmp)
  dev.off()
}







