### full pipeline 


#1. run transcriptomics
## change here to run once for mirs and once for rnas 
source(paste0('ppmi/setup_os.R'))

## SET VISIT AND PROCESS_MIRNAS here 
source(paste0(script_dir,'ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))
source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))
source(paste0(script_dir,'ppmi/deseq_analysis.R'))


# SET VISIT HERE FOR THE SCRIPT 
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir,'ppmi/olink-preprocessing_ppmi_se.R'))
source(paste0(script_dir,'ppmi/DE_tutorial_ppmi.R'))


