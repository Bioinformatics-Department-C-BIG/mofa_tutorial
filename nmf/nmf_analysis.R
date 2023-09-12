





#### Correlations ####


library(psych)


dim(covariates)
cor <- psych::corr.test(covariates,h_t, method = "pearson", adjust = "BH")


plot_covars<-function(cor, covars_to_plot ,fname, plot='r', labels_col=FALSE){
  if (is.logical(labels_col)){labels_col=covars_to_plot}
  
  
  labels_col<-ll$V2[match(covars_to_plot, ll$V1)]
  
  if (plot=="r") {
    
    stat <- cor$r
    png(paste0( outdir_nmf, '/',fname,  '.png' ))
    corrplot::corrplot(stat[covars_to_plot,], tl.col = "black", title="Pearson correlation coefficient")
    dev.off()  
    
    write.csv(stat[covars_to_plot,], paste0(outdir_nmf, 'cor_',mod, '.csv'  ))
    
    
  } else if (plot=="log_pval") {
    stat <- cor$p
    #stat[stat>alpha] <- 1.0
    if (all(stat==1.0)) stop("All p-values are 1.0, cannot plot the histogram")
    stat <- -log10(stat)
    
    
    stat_filt<-  as.data.frame(stat[covars_to_plot,])
    stat_filt<-as.data.frame(na.omit(stat_filt))
    stat[is.infinite(stat)] <- 1000
    #if (transpose) stat <- t(stat)
    #if (return_data) return(stat)
    col <- colorRampPalette(c("lightgrey", "red"))(n=100)
    
    png(paste0( outdir_nmf, '/', fname, '.png' ), res=300, width=10, height=7, units='in')
    
    pheatmap::pheatmap(t(stat_filt), main="log10 adjusted p-values", cluster_rows = TRUE, color=col, 
                       labels_col = labels_col )
    
    
    dev.off()  
  }
}




colnames(cor$r) 
rownames(cor$sef)
cors_non_na<-names(which(!is.na( cor$p[,1])))
T<--log10(0.00)


to_remove_regex<-'DATE|REC_ID|UPDATE|ORIG_ENTR|INFO|PATNO|cluster'
to_remove_covars<-grepl( to_remove_regex, rownames(cor$p))


sig<-which( rowMins(cor$p)<0.01 & rowMaxs(cor$r)<0.95 & !to_remove_covars ) 
sig<-which( rowMins(cor$p)<0.01 & !to_remove_covars ) 


sig2<-which( rownames(cor$p) %in% selected_covars2 ) 
rownames(cor$p)[which( rownames(cor$p) %in% selected_covars2 ) ]


sig_broad<-which( rownames(cor$p) %in% selected_covars_broad &  !(rownames(cor$p) %in% 'con_putamen')  ) 


plot='r'



fname='covariates'
covars_to_plot=sig
sig

plot_covars(cor, sig, fname='covariates_sig', plot='r')

plot_covars(cor, sig2, fname='covariates',  plot='r')
plot_covars(cor, sig_broad, fname='covariates_broad', plot='r')


plot_covars(cor, sig, fname='covariates_sig_logpval', plot='log_pval')
plot_covars(cor, sig2, fname='covariates_logpval',  plot='log_pval')
plot_covars(cor, sig_broad, fname='covariates_broad_logpval',  plot='log_pval')

  
#} else {
#  stop("'plot' argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
#}

which( rowMins(cor$p)<0.05)

cor$r['CONCOHORT',]
cor$p['CONCOHORT',]

sel_factors_s<-which(cor$p['COHORT',]<0.05)

stat[sig,]

#### Top Weights ####


ll
w<-basis(res)
write.csv(w, paste0(nmf_outdir, '/top_weights/weights_all_factors.csv'))

for (fn in 1:NFACTORS){
  w_ordered<-w[, fn][order(abs(w[, fn]), decreasing = TRUE)]
  write.csv(w_ordered, paste0(nmf_outdir, '/top_weights/weights', fn, '.csv'))
  
  
}







#### Cluster ####
# 1. K-means 
#  TODO: other clusters 
#sel_factors=c(2,4)
#sel_factors=c(1,2,3)
k_centers=3  ## best one 
k_centers=6  ## best one 

clusters_single <- kmeans(t(h)[,sel_factors_s], centers = k_centers)
#clusters_single_rna<-clusters_single

covariates$cluster_s<-clusters_single$cluster[match(rownames(covariates),names(clusters_single$cluster))]
if (mod=='RNA'){covariates$cluster_s_rna<-clusters_single}
covariates$cluster_m<-clusters_mofa$cluster[match(rownames(covariates),names(clusters_mofa$cluster))]
covariates$clusters_mofa_moca<-clusters_mofa_moca$cluster[match(rownames(covariates),names(clusters_mofa_moca$cluster))]

chisq.test(df1$cluster_m, df1$COHORT)


clusters_single$size
clusters_mofa$size

chisq.test(clusters_single$cluster,covariates$COHORT )
chisq.test(clusters_single$cluster,covariates$COHORT )

MutInf(clusters_single$cluster,covariates$COHORT)



df1=covariates
chisq.test(df1$cluster_m, df1$COHORT)

chisq.test(df1$cluster_s, df1$COHORT)
to_test<-selected_covars2
selected_covars2
to_test<-selected_covars2[!selected_covars2 =='PDSTATE']

to_test<-selected_covars_broad[!selected_covars_broad %in% c('STAIAD26')]
to_test=c(to_test, 'td_pigd', 'td_pigd_on', 'td_pigd_old', 'COHORT')


get_pval_cor<-function(all_tests){
  
  
  as.data.frame(lapply(all_tests, function(x){
    return(x$p.value)
  }))}


check_if_in_metadata<-function(to_test, df1){
  #' check if inside the colnames and not all is NA 
  #'  @param df1 df to test description
  #' @param to_test description
  all_na<-names(which(apply(df1, 2, function(x){all(is.na(x))})))
  to_test=to_test[to_test %in% colnames(df1) & !(to_test%in% all_na)]
  return(to_test)
}

get_cluster_assoc<-function(df1, cluster_col, to_test){
  
  
  to_test<-check_if_in_metadata( to_test,df1)
  all_tests<-apply(df1[,to_test],2,  kruskal.test, g=as.factor(df1[,cluster_col ] ), na.rm=TRUE)
  all_tests_pval<-get_pval_cor(all_tests)
  return(all_tests_pval)
}

cluster_col='cluster_s'


get_cluster_mi<-function(df1, cluster_col,to_test){
  
  to_test<-check_if_in_metadata( to_test,df1)
  
  all_tests<-apply(df1[,to_test],2,  MutInf, y=as.factor(df1[,cluster_col] ))
  return(all_tests)
}

all_tests_m<-get_cluster_assoc(df1, 'cluster_m',to_test) 
all_tests_s<-get_cluster_assoc(df1, 'cluster_s',to_test) 
all_tests_moca<-get_cluster_assoc(df1, 'clusters_mofa_moca', to_test) 

all_tests_m_mi<-get_cluster_mi(df1, 'cluster_m', to_test) 
all_tests_s_mi<-get_cluster_mi(df1, 'cluster_s', to_test) 


all_tests_compare<-rbind(all_tests_m, all_tests_s);rownames(all_tests_compare)<-c('MOFA', 'RNA')
#all_tests_compare<-rbind(all_tests_compare, all_tests_moca); rownames(all_tests_compare)<-c('MOFA', 'RNA', 'mofa_moca')


### Associations with Clusters all_tests_compare
labels_col<-mt_kv$V2[match(colnames(all_tests_compare),mt_kv$V1)]

#### Evaluation #### 
png(paste0(output_files, 'cluster_covars.png'), height=5, width=8, units='in', res=600)
ph<-pheatmap(as.matrix(-log10(all_tests_compare)), cluster_rows = FALSE, labels_col=labels_col)
#pheatmap(as.matrix(all_tests_compare), cluster_rows = FALSE)
ph
dev.off()


all_tests$COHORT$p.value
all_tests$COHORT$p.value

kruskal.test(df1$NHY, as.factor(df1$cluster_m ))
kruskal.test(df1$NHY, as.factor(df1$cluster_s ))

kruskal.test(df1$NP3_TOT, as.factor(df1$cluster_m ))

kruskal.test(df1$NP2_TOT, as.factor(df1$cluster_m ))
kruskal.test(df1$NP2_TOT, as.factor(df1$cluster_s ))




dim(df1)







