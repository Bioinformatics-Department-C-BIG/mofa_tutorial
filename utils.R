


library('plyr')

library('dplyr')
group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*neural.*|.*deep.*|.*autoencoder.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                   '.*pca.*|.*cluster.*|.*kmeans.*|.*pins.*|.*movis.*|.*k means.*|.*partitioning.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*', 
                   '.*snf.*|.*net.*|.*piumet.*|.*omics integrator.*',
                   '.*gsea.*|.*enrichment.*', 
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*kernel.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*',
                   '.*david.*|.*metaboanalyst.*|.*nemo.*|.*adas.*'), 
                 c( "ML/DL Classification", 
                    'ML/DL clustering',
                    'regression', 
                    'factor analysis', 
                    'network', 
                    'enrichment',
                    'correlation-based',
                    'correlation',
                    'kernel learning', 
                    'partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}



group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*neural.*|.*deep.*|.*autoencoder.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                   '.*pca.*|.*cluster.*|.*kmeans.*|.*pins.*|.*movis.*|.*k means.*|.*partitioning.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*', 
                   '.*snf.*|.*net.*|.*piumet.*|.*omics integrator.*',
                   '.*gsea.*|.*enrichment.*', 
                   '.*cca.*|.*smccnet.*|.*canonical correlation analysis.*',
                   '.*correlation.*', 
                   '.*kernel.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*',
                   '.*david.*|.*metaboanalyst.*|.*nemo.*|.*adas.*'), 
                 c( "ML/DL Classification", 
                    'ML/DL clustering',
                    'regression', 
                    'factor analysis', 
                    'network', 
                    'enrichment',
                    'correlation-based',
                    'correlation',
                    'kernel learning', 
                    'partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}


group_methods<-function(df, Var){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                   ".*autoencoder.*|.*umap.*|.*tsne.*|.*deep.*.|*neural.*|.*graph convolutional network.*|.*cdrscan.*|.*deepomix.*",
                   '.*cluster.*|.*kmeans.*|.*pins.*|.*k means.*|.*kernel.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*|.*partitioning.*|.*pca.*|.*diverse.*|.*pathme.*', 
                   '.*net.*|.*piumet.*|.*omics integrator.*|.*inet.*',
                   '.*snf.*|.*coni.*|.*netdx.*|.*nem-tar.*|.*paradigm.*',
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*|.*panther.*|.*david.*|.*gsea.*|.*enrichment.*',
                   '.*metaboanalyst.*|.*nemo.*|.*adas.*|.*movics.*|.*mousse.*|.*timeg.*|.*miodin.*|.*ioda.*'), 
                 c( "ML Classification", 
                    'JDR - NonLinear',
                    'JDR - Linear',
                    'Regression', 
                    'JDR - Linear - MF', 
                    'Network analysis',
                    'Graph-based', 
                    'JDR - correlation-based',
                    'correlation',
                    'JDR - Linear - PLS',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}




relabel_objectives_short<-function(df_to_plot){
  df_to_plot$labels<-df_to_plot$key_names
  ind<-df_to_plot$labels%in% c('connect molecular patterns to phenotypic traits')
  df_to_plot[ind,]$labels<-'connect molecular patterns to \n phenotypic traits'
  df_to_plot[ind,]$key_names<-'connect molecular patterns to \n phenotypic traits'
  
  ind<-df_to_plot$labels%in% c('understand molecular mechanisms')
  df_to_plot[ind,]$labels<-'understand \n molecular mechanisms'
  df_to_plot[ind,]$key_names<-'understand \n molecular mechanisms'
  
  
  return(df_to_plot)
}

relabel_objectives<-function(obj_col){
  
  obj_col <- mapvalues(obj_col, from=c("understand molecular mechanisms",
                                       "understand regulatory processes",   
                                 "connect molecular patterns to phenotypic traits",
                                 "Diagnosis/Prognosis",
                                 "biomarker discovery",
                                 
                                 "subtype identification",
                                 "multiomics pathway analysis",
                                 
                                 "drug response prediction"), 
                     to=c("UN","UN", "CO", "DI", "BM", "SI","PA" ,  "DR"))
    return(obj_col)
}




group_objectives_method<-function(df, Var1){
  #'Group objective code column 
  #'These groups are for objective - method
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand.*'),
                                  c('Diagnosis/Prognosis', 'understand molecular mechanisms')))
  return(df)
}


group_disease<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                c(".*alzheimer.*|.*amyotrophic.*|.*anxiety*|.*depression.*|.*parkinson.*|.*autism.*",
                  '.*cardio.*|.*heart.*|.*coronary.*',
                  '.*bowel.*|.*hep.*|.*liver.*|.*nafld.*|.*crohn.*', 
                  '.*arthritis.*|.*osteo.*', 
                  '.*diabetes.*', 
                  '.*pulmonar.*|.*lung.*|.*copd.*|.*smoking.*',
                  '.*cancer.*|.*carcinoma'), 
                c('Nervous system', 
                  'Cardiovascular',
                  'Gastrointestinal', 
                  'Musculoskeletal', 
                  'Endrocrine', 
                  'Pulmonary', 
                  'Cancer'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}
