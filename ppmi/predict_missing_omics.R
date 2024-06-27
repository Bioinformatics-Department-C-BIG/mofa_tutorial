


  #  mofa imputation 
# remove for one sample one modality
# impute it 
# and test the corelation 

    
  

library(ggplot2)

# TODO: check again if the DE genes/ proteins are also well predicted!!
# Example on the CLL data
filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
MOFA_CLL <- MOFA2::load_model(filepath)

# predict drug response data using all factors
view="proteomics_t_plasma"
predicted_proteins <- MOFA2::predict(MOFAobject, view=,groups=1,add_intercept = FALSE)
#predicted_proteins <- MOFA2::predict(MOFAobject, view="proteomics_csf",groups=1,add_intercept = FALSE)

predicted_proteins

pr_patients_i = 2
mofa_multi_to_use_train<-mofa_multi

library(MultiAssayExperiment)

# Identify the specific assay and sample
assay_name <- "proteomics_csf"
sample_name <- 1





# added se filter   
# 
view = 'proteomics_plasma'
data_full[[view]]


sample_name
data_full_train<-data_full
data_full_train[[view]][,pr_patients_i]<-NA
data_full[[view]][,pr_patients_i]

colnames(data_full_train[[view]])

all_prot<-data_full_train[[view]] 

TOP_N =0.1

hvp<-selectMostVariable(all_prot, TOP_N)
high_var_prot<-rownames(highly_variable_proteins_predict)

sid<-'12224_V08'

sid
mofa_multi_to_use_train<-create_multi_experiment(data_full_train, combined_bl_log)

prot_samples<-colnames(assays(mofa_multi_to_use_train)[[view]])

# subset the test sample
prot_patient_samples<-mofa_multi_to_use_train[,mofa_multi_to_use_train$COHORT%in% 1]$PATNO_EVENT_ID
prot_samples<-intersect(prot_patient_samples,prot_samples)
s_ind<-10
sid = prot_samples[s_ind]
all_predictions = list()
for (sid in prot_samples){




      # remove from data_full and recreate multi assay 

      original_test<-data_full[[view]][,sid]

      data_full_train<-data_full
      data_full_train[[view]][,sid]<-NA # remove proteomics for the test sample
      mofa_multi_to_use_train<-create_multi_experiment(data_full_train, combined_bl_log)

      mofa_multi_to_use_train






      outdir_train = paste0(outdir, '../', sid)

      outdir_train = paste0(outdir_orig,shorten_path(outdir_full),'_',view,sid )



      MOFAobject_train=run_mofa_get_cors(mofa_multi_to_use_train, N_FACTORS, force=FALSE, outdir=outdir_train)





      #### 
      # Leave one out cross validation:
      # Now predict on hold-out 

      sm<-samples_metadata(MOFAobject_train)

      # data.frame(get_data(MOFAobject, views=view, groups = 1)
      # all predictions for selected view
      predicted_all<- MOFA2::predict(MOFAobject_train, view=view)[[1]]$group1
      # prediction for selected sample 
      predicted_test<-as.data.frame(predicted_all)[,sid]




      colnames(as.data.frame(get_data(MOFAobject, views=view)[[1]]$group1))

      length(original_test)
      original_test[high_var_prot]
      predicted_test[high_var_prot]
      names(predicted_test)<-names(original_test)
      names(predicted_test)

      cor_predicted_original<-corr.test(original_test[high_var_prot],predicted_test[high_var_prot] )


      all_predictions[[sid]]<-c(cor_predicted_original$r, cor_predicted_original$p.adj)
}

#all_predictions<-as.data.frame(all_predictions)
write.csv(all_predictions, paste0(outdir,'/../all_predictions.csv'))






original_all<-rownames(na.omit(t(MOFAobject@data[[view]][[1]])))
original_all<-t(na.omit(t(MOFAobject@data[[view]][[1]])))

non_na<-colnames(original_all)
s_id='4029_V08'
predicted=predicted_proteins[[ view]][[1]][, non_na]
original<-MOFAobject@data[[ view]][[1]][,non_na]
original
sm<-samples_metadata(MOFAobject)


sm_coh<-sm[sm$PATNO_EVENT_ID%in% non_na,]$COHORT
sm[sm_coh==1,'PATNO_EVENT_ID']

non_na

original_ps<-original[sm_coh%in% c(2),]
predicted_ps<-predicted[sm_coh%in% c(2),]

# corelations after selection of HC, or PD separately
cors_train<-sapply(seq.int(dim(original_ps)[2]), function(i) cor(original_ps[,i], predicted_ps[,i]))
original_ps
mean(sapply(seq.int(dim(original_ps)[2]), function(i) cor(original_ps[,i], predicted_ps[,i])))
hist(cors_train)

corr.test(original, predicted)
# what about on the DE genes?
# predict all views using all factors (default)
predictedAll <- MOFA2::predict(MOFAobject)

predictedAll
# predict Mutation data using all factors, returning Bernoulli probabilities
predictedMutations <- predict(MOFA_CLL, view="Mutations", type="response")

# predict Mutation data using all factors, returning binary classes
predictedMutationsBinary <- predict(MOFA_CLL, view="Mutations", type="inRange")

# Compare the predictions with the true data
pred <- as.numeric(predictedAll$Drugs)
true <- as.numeric(getTrainData(MOFA_CLL)$Drugs)
qplot(pred,true) + geom_hex(bins=100) + coord_equal() + 
   geom_abline(intercept=0, slope=1, col="red")

# Example on the scMT data
filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
MOFA_scMT <- loadModel(filepath)

# Predict all views using all factors (default)
predictedAll <- predict(MOFA_scMT)
 
# Compare the predictions with the true data
view <- "RNA expression"
pred <- as.numeric(predictedAll[[view]])
true <- as.numeric(getTrainData(MOFA_scMT)[[view]])
qplot(pred,true) + geom_hex(bins=100) + coord_equal() + 
   geom_abline(intercept=0, slope=1, col="red") 



