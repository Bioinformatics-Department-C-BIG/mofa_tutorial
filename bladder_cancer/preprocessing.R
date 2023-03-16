
sysinf <- Sys.info()


####################################
#### Preliminary analysis with PCA
os <- sysinf['sysname']
getwd()
if ( os  == 'Darwin'){
  dir='bladder_cancer/data/'
}else if ( os ==   'Windows'){
  dir='bladder_cancer/data/'
}else{
  dir='/data8TB/efiath/git/mofa/'
}
output='bladder_cancer/plots/'

transpose_matrix<- function(df.aree){
  n <- df.aree$Symbol
  # transpose all but the first column (name)
  df.aree <- as.data.frame(t(df.aree[,-1]))
  colnames(df.aree) <- n
  #df.aree$myfactor <- factor(row.names(df.aree))
  return(df.aree)
  
}



preprocess_raw_data<-function(df, cut_n){
  
  
  # first remove all zero entries
  if (cut_n){ df<-df[,1:cut_n] } 
  
  df<- as.data.frame(apply(df, 2, function(x) as.numeric(x)) )
  ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
  df<-df[,!ind]
  
  return(df)
  
  
  
}

selectMostVariable<-function(vsn_mat,q){
  #' selects rows ie. genes must be rows
  #' Selects top q most variable genes
  #' Ideally take vsn transformed dataset!
  #' @param: vsn_mat: genes/proteomics matrix after vsn/vst transform
  #' q: top q genes/
  variances <- apply(vsn_mat, 1, var, na.rm=TRUE)
  topx<-names(variances[order(variances, decreasing = TRUE)])[1:round(length(variances)*q, digits=0)]
  vsn_mat <- vsn_mat[topx, ]
  NROW(vsn_mat);dim(vsn_mat)
  if (is.null(topx)){print('Warning: zero most variable features returned')}
  return(vsn_mat)

  }


filter_most_var<-function(df,most_var, ng){  # take the most variable entries 
  if (most_var){
    n=round(dim(df)[2]/ng)
    mads<-apply(df,2,mad)
    df_selected=df[,rev(order(mads))[1:n]]
    return(df_selected)
  }else {
    return(df)
  }
}



X1_raw<-read.csv(file = paste0(dir,'RNAseq_BladderCancer.csv' ))
X2_raw<-read.csv(file = paste0(dir,'Proteomics_BladderCancer.csv' ))
Y_raw<-read.csv(file = paste0(dir,'pheno_BladderCancer.csv' ), nrows = 16)


X1_t_raw<-transpose_matrix(X1_raw)
X2_t_raw<-transpose_matrix(X2_raw)


most_var=TRUE
X1_t_cut<-preprocess_raw_data(X1_t_raw, cut_n=27000)
X2_t_cut<-preprocess_raw_data(X2_t_raw, cut_n = FALSE)



ng_p=70
ng_g=25


q=0.7
##### Select most variable genes
