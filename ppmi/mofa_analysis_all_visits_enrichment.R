
### continuation from mofa analysis all visits 
#source('')
## problem: how can I know which factors are relevant ?? 
## since I ma not doing a controls?? 

f1<-all_fs_merged1[str_detect(all_fs_merged1[,2], 'PARKINSON'),'factor']
f2<-all_fs_merged2[str_detect(all_fs_merged2[,2], 'PARKINSON'),'factor']

all_f_park<-c(f1,f2)









####
###### SAVE ONLY HIGH VAR RNA FOR CLUEGO
### that are also related to parkinsons?? 

variances<-get_variance_explained(MOFAobject, views =c(3) )$r2_per_factor$V06
high_var_fs<-which(variances[,'RNA']>1)
high_var_fs
view='RNA'

all_f_park
cluego2<-paste0(outdir, 'top_weights_vals_by_view_CLUEGO_park_factor_', view, '_T_', T, '.txt')

#### ONLY SAVE HIGH VARIANCE in RNAS
all_weights2<-MOFA2::get_weights(MOFAobject,
                                 views = view,
                                 factors=all_f_park,
                                 as.data.frame =TRUE)  
# threshold each? 

all_weights_filt<-all_weights2[abs(all_weights2$value)>T,]
ens_ids<-gsub('\\..*', '', all_weights_filt$feature)
write.csv(ens_ids,cluego2,
          row.names = FALSE, quote=FALSE)


