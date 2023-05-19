
#script_dir='/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/'
#os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))
source(paste0(script_dir, '/setup_os.R'))

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
#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)


#### TODO: SAVE SE FILT SO WE CAN RELOAD in the next script 

if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
#source(paste0(script_dir,'/../bladder_cancer/preprocessing.R'))
source(paste0(script_dir,'/utils.R'))

output_1=paste0(data_dir,'ppmi/plots/proteomics/')
outdir_orig<-paste0(data_dir,'ppmi/plots/')
output_files<-paste0(data_dir,'ppmi/output/')


metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)
combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)





metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)


### TODO: filter out the cohort too before processig !! 

VISIT=c('BL')
process_mirnas=FALSE
source(paste0(script_dir, '/config.R' ))

param_str<-paste0(TOP_PN)


run_vsn
NORMALIZED




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

outdir<-outdir_orig

#### Read in 
prot_bl_wide_unlog<-as.matrix(fread(in_file_original, header=TRUE), rownames=1)
in_file_original
proteomics<-prot_bl_wide_unlog
  
#proteomics[2,]


colnames(proteomics)


#### FILTERING LOW VALUES 
# Remove rows with 90% NA 

#### TODO: move this after the selection of the cohort!!! 
df<-proteomics
proteomics <- df[rowSums(is.na(df)) < round(0.2*ncol(df)), ]
dim(df); dim(proteomics)


### filter here before editing more 
## filter out rows with very low min count
df<-proteomics; dim(df)
min.count= quantile(df, na.rm = TRUE, 0.01)
min.count= min(df, na.rm = TRUE)
min.count
dim(df)
### KEEP
 ## kEEP THE ROWS THAT HAVE MORE THAN 80% NON NA VALUES 
keep <- rowSums(df>min.count, na.rm = TRUE) >= round(NA_PERCENT*ncol(df))
length(which(keep))
proteomics<-proteomics[keep,]
dim(df); dim(proteomics)

##keep <- colSums(df>min.count, na.rm = TRUE) >= round(0.9*ncol(df))
#length(keep)
#length(which(!keep))
#proteomics<-proteomics[,keep]




#### MAKE NUMERIC 
raw_counts_all=proteomics
class(raw_counts_all) <- "numeric"
## They seem to have taken averages for replicas so need to fix 
#raw_counts_all<-round(raw_counts_all)
dim(proteomics)
dim(raw_counts_all)
head(raw_counts_all[,2])

#### Input to se 




dim(proteomics)
data<-proteomics
data<-as.data.frame(data)

data$name<-c(rownames(data))
data$ID<-data$name
data_columns=seq(1:dim(proteomics)[2])

#sample<-colnames(proteomics_se)
head(sample)




#exp_design = data.frame(label=sample,
#                        condition=sample, 
#                        replicate=rep(1, dim(proteomics)[2]))
#exp_design

#se <- make_se(data, data_columns,exp_design)
## just a quick filter because it takes too much memory
#raw_counts_all_by_visit<-raw_counts_all[,grep(VISIT, colnames(raw_counts_all ))]
proteomics_se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined)
dim(proteomics_se)

head(is.nan(as.matrix(data)))
boxplot(log(data[1:16]))
head(raw_counts_all[2,])

#View(head(assays(proteomics_se)[[1]]))


cbind(proteomics_se$COHORT_DEFINITION,proteomics_se$COHORT  )

sel_coh
##### filter here by visits
se_filt<-proteomics_se[,(proteomics_se$EVENT_ID %in% VISIT & proteomics_se$COHORT %in% sel_coh )]

se_filt<- filter_se(proteomics_se, VISIT, sel_coh, sel_ps)

dim(se_filt)
table(se_filt$COHORT)

cbind(se_filt$COHORT_DEFINITION,se_filt$COHORT  )

### TODO: save se filt here : with or without VISIT included..? 
Sample<-colnames(se_filt)
sample_info<-DataFrame(Sample=Sample)


tmp<- assays(se_filt)[[1]]
dim(tmp)

# Filter and normalize
### ERROR: Error in vsnML(sv) : L-BFGS-B needs finite values of 'fn'
#normalized_data1<-assays(normalize_vsn(se_filt))[[1]]
#normalized_data1<-normalize_vsn(se_filt)
#meanSdPlot(normalized_data1)
#vsn_mat<-assays(normalized_data1)[[1]]

#ggsave(paste0(output_1,'meansd_normalize_vsn_', p_params_out,'.png' ))

#normalized_data<-varianceStabilizingTransformation(tmp  )
### TODO: filter before normalization!!! 
normalized_data<-justvsn(tmp)
vsn::meanSdPlot(normalized_data)



ggsave(paste0(output_1,'meansd_justvsn_', p_params_out,'.png' ), width = 5, height=3)
# Check plot after vsn
#View(normalized_data)


vsn_mat<-normalized_data

# Select the top most variable proteins
## TODO: fix the bug in selectMostVariable
TOP_PN
for (most_var in c(0.05, 0.1,0.15,0.2,0.25,0.3,  0.9,0.75,0.5)){
     # TODO: IF YOU CHANGE THIS MAKE SURE TO CHANGE THE ONE IN MOFA 
    #p_params<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s, sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
  
    p_params_out<- paste0(VISIT_S, '_',TISSUE, '_', most_var, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
    highly_variable_proteins_outfile<-paste0(output_files, p_params_out , '_highly_variable_proteins_mofa.csv')
  
  
    highly_variable_proteins_mofa=selectMostVariable(vsn_mat, most_var)
    
    write.csv(highly_variable_proteins_mofa,highly_variable_proteins_outfile)
    print(dim(highly_variable_proteins_mofa))
}
png(paste0(output_1,'hist_high_var_', p_params_out,'.png' ))
hist(highly_variable_proteins_mofa)
dev.off()


png(paste0(output_1,'hist_', p_params_out,'.png' ))
hist(vsn_mat)
dev.off()



#dev.off()
colnames(highly_variable_proteins_mofa)




#### save all and load 
##deseq2Results <- results(deseq2Data, contrast=c('COHORT', 1,2))
datalist=list( vsn_mat, se_filt)
saveRDS(datalist,prot_vsn_se_filt_file)



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


highly_variable_proteins_outfile

