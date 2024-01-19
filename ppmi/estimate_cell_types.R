


library(granulator)
# load datasets for deconvolution of PBMC RNA-seq data
load_ABIS()


# print TPM values in bulk RNA-seq
bulkRNAseq_ABIS[1:5, 1:5]

# print TPM values in reference profile matrix
sigMatrix_ABIS_S0[1:5, 1:5]


# print measured cell type proportions (percentages)
groundTruth_ABIS[1:5, 1:5]



# create list if multiple signature matrices to test simultaneously
sigList = list(
  ABIS_S0 = sigMatrix_ABIS_S0,
  ABIS_S1 = sigMatrix_ABIS_S1, 
  ABIS_S2 = sigMatrix_ABIS_S2, 
  ABIS_S3 = sigMatrix_ABIS_S3)



  # plot signature matrix similarity matrices
plot_similarity(sigMatrix=sigList)

vst_cor_all_vis_data<-loadRDS(vst_cor_all_vis)


process_mirnas=FALSE; 
source(paste0(script_dir, '/ppmi/config.R'))

se_rnas<-load_se_all_visits(input_file = input_file, combined=combined_bl_log); 

bulkRNAseq<-assay(se_rnas)
rownames(bulkRNAseq)
rownames(bulkRNAseq)<-get_symbols_vector(rownames(bulkRNAseq))


# deconvolute input data using all available methods by default
decon <- deconvolute(m = bulkRNAseq, sigMatrix = sigList, methods=get_decon_methods()[1:4])


estimated<-decon$proportions$nnls_ABIS_S0
max(decon$proportions$nnls_ABIS_S0)

# plot cell type proportions for svr model on ABIS_S0 reference profile
plot_proportions(deconvoluted = decon, method = 'nnls', signature = 'ABIS_S2')
# plot cell type proportions
plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)

ground_truth<-as.data.frame(colData(se_rnas))

ground_truth$PATNO_EVENT_ID
match(rownames(estimated),ground_truth$PATNO_EVENT_ID )


match(ground_truth$PATNO_EVENT_ID , rownames(estimated))


estim_matched<-estimated[match(ground_truth$PATNO_EVENT_ID , rownames(estimated)),]

estim_matched$Neutrophils.LD

sm<-samples_metadata(MOFAobject)
se_rnas$Neutrophil.Score

library('psych')
corr.test(cbind(estim_matched$Neutrophils.LD, ground_truth$`Neutrophils....`))

corr.test(estim_matched$Neutrophils.LD, ground_truth$`Neutrophils....`)$p
corr.test(estim_matched$Neutrophils.LD, ground_truth$`Neutrophils....`)$p
corr.test(estim_matched$Neutrophils.LD, ground_truth$Neutrophil.Score)$p
estimates<-