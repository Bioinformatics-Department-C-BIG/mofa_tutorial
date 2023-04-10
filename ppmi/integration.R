
#### 
#Calculate the Pearson correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.

BiocManager::install('miRLAB')
library('miRLAB')


library(miRLAB)


dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
results=Pearson(dataset, 1:3, 4:18) 

mirs_vals<-assays(mofa_multi)[[1]]
l1<-dim(mirs_vals)[2]
rnas_vals<-assays(mofa_multi)[[2]]

dataset1<-rbind( mirs_vals, rnas_vals)


d2<-t(dataset1)
d2_f<-paste0(out_compare,'file1.csv')
write.csv(d2,d2_f, row.names = FALSE)

results=Pearson(d2_f, dim(mirs_vals)[1],dim(rnas_vals)[1]) 

outdir_s_1<-paste0(outdir_orig, '/single/', param_str_g_f, des)


