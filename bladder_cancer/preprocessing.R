
sysinf <- Sys.info()


####################################
#### Preliminary analysis with PCA
os <- sysinf['sysname']
if ( os  == 'Darwin'){
  dir='/Users/efiathieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'
}else if ( os ==   'Windows'){
  dir='E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'
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




ng_p=round(3/2,2)
ng_g=round(10,2)
X1_t_cut<-preprocess_raw_data(X1_t_raw, cut_n=27000)
X1_t_most_var<-filter_most_var(X1_t_cut,most_var,ng_g)


X2_t_cut<-preprocess_raw_data(X2_t_raw, cut_n = FALSE)
X2_t_most_var<-filter_most_var(X2_t_cut, most_var,ng_p)


X1_t_1<-as.data.frame(cpm(X1_t_most_var, log = TRUE) )
X2_t_1<-as.data.frame(cpm(X2_t_most_var, log = TRUE) )

