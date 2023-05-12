
#### 
#Calculate the Pearson correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.

#BiocManager::install('miRLAB')
library('miRLAB')
library(miRLAB)
out_compare<-'ppmi/plots/single/compare/'

mofa_multi_complete
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
results=Pearson(dataset, 1:3, 4:18) 

mirs_vals<-assays(mofa_multi_rna_mir_complete)[[1]]
l1<-dim(mirs_vals)[2]
dim(mirs_vals)
l1
rnas_vals<-assays(mofa_multi_rna_mir_complete)[[2]]
dim(rnas_vals)
dataset1<-rbind( mirs_vals, rnas_vals)
dim(dataset1)

d2<-t(dataset1)
d2_f<-paste0(out_compare,'file1_', TOP_GN, '_', TOP_MN,'.csv')
write.csv(d2,d2_f, row.names = FALSE)

results=Pearson(d2_f, 1:dim(mirs_vals)[1],(dim(mirs_vals)[1]+1):dim(dataset1)[1])


hist(results)
anticor= results <= -0.35

write.csv(results, paste0(outdir_s, '/', TOP_GN, '_', TOP_MN, '_', 'cor_results.csv') )
#outdir_s_1<-paste0(outdir_orig, '/single/', param_str_g_f, des)

TOP_GN
TOP_MN
