
#script_dir='/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/'
#os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'
getwd()
source('ppmi/setup_os.R')


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
source(paste0(script_dir,'ppmi/utils.R'))
output_1=paste0(data_dir,'ppmi/plots/proteomics/')
outdir_orig<-paste0(data_dir,'ppmi/plots/')
output_files<-paste0(data_dir,'ppmi/output/')





combined_bl_log<-load_metadata()


### TODO: filter out the cohort too before processig !! 

all_visits = c('BL', 'V02', 'V04', 'V06', 'V08')

for (VISIT in all_visits){
    #process_mirnas=FALSE
    #TISSUE='CSF'

    source(paste0(script_dir, 'ppmi/config.R' ))


    param_str<-paste0(TOP_PN)
    p_params_in
    pr_ext<-'.csv'
    if (NORMALIZED){
      pr_ext<-  '_no_log.csv'
      # if we dont run vsn we want to take the log values 
      if (!run_vsn){
        pr_ext<-'.csv'
      }
    }else{
        pr_ext<-  '_no_log.csv'
    }
    in_file_original<-paste0(output_files, 'proteomics_', p_params_in, pr_ext)


    #### Read in 
    prot_bl_wide_unlog<-as.matrix(fread(in_file_original, header=TRUE), rownames=1)

    proteomics<-prot_bl_wide_unlog
    #### FILTERING LOW VALUES 
    # Remove rows with 90% NA 
    #### TODO: move this after the selection of the cohort!!! 

    raw_counts_all<-pre_process_proteomics(proteomics)
    proteomics_se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined_bl_log)
    dim(proteomics_se)

    ##### filter here by visits and cohort
    se_filt_proteins<- filter_se(proteomics_se, VISIT, sel_coh, sel_ps)

    ### TODO: save se filt here : with or without VISIT included..? 
    Sample<-colnames(se_filt_proteins)
    sample_info<-DataFrame(Sample=Sample)

    tmp<- assays(se_filt_proteins)[[1]]


    # Filter and normalize
    ### ERROR: Error in vsnML(sv) : L-BFGS-B needs finite values of 'fn'

    ### TODO: filter before normalization!!! 


    # Select the top most variable proteins
    ## TODO: fix the bug in selectMostVariable

    p_params_out<- paste0(VISIT_S, '_',pr_project_id, '_', TISSUE, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
    p_params_out
    normalized_data<-justvsn(tmp)
    vsn::meanSdPlot(normalized_data)



    ggsave(paste0(output_1,'meansd_justvsn_', p_params_out,'.png' ), width = 5, height=3)


    vsn_mat<-normalized_data
    head(vsn_mat)

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


}






