





#### Correlations ####


library(psych)


dim(covariates)
cor <- psych::corr.test(covariates,h_t, method = "pearson", adjust = "BH")



colnames(cor$r) 
rownames(cor$sef)
cors_non_na<-names(which(!is.na( cor$p[,1])))
T<--log10(0.00)


to_remove_regex<-'DATE|REC_ID|UPDATE|ORIG_ENTR|INFO|PATNO'
to_remove_covars<-grepl( to_remove_regex, rownames(cor$p))


sig<-which( rowMins(cor$p)<0.01 & rowMaxs(cor$r)<0.95 & !to_remove_covars ) 
sig<-which( rowMins(cor$p)<0.01 & !to_remove_covars ) 


sig<-which( rownames(cor$p) %in% selected_covars2 ) 
sig





which( rowMins(cor$p)<0.05)

cor$r['CONCOHORT',]
stat[sig,]


plot='r'


chosen_covars
  
if (plot=="r") {
  stat <- cor$r
  png(paste0( outdir_nmf, './covariates.png' ))
  
  corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")
  dev.off()  
  
  write.csv(stat[sig,], paste0(outdir_nmf, 'cor_',mod, '.csv'  ))
} else if (plot=="log_pval") {
  stat <- cor$p
  stat[stat>alpha] <- 1.0
  if (all(stat==1.0)) stop("All p-values are 1.0, cannot plot the histogram")
  stat <- -log10(stat)
  
  
  stat_filt<-  as.data.frame(stat[sig,])
  stat_filt<-as.data.frame(na.omit(stat_filt))
  stat[is.infinite(stat)] <- 1000
  if (transpose) stat <- t(stat)
  if (return_data) return(stat)
  col <- colorRampPalette(c("lightgrey", "red"))(n=100)
  
  png(paste0( outdir_nmf, './covariates_logpval.png' ), res=300, width=10, height=7, units='in')
  
  pheatmap::pheatmap(t(stat_filt), main="log10 adjusted p-values", cluster_rows = TRUE, color=col)
  
  dev.off()  
  
  
} else {
  stop("'plot' argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
}


#### Top Weights ####



w<-basis(res)
write.csv(w, paste0(nmf_outdir, '/top_weights/weights_all_factors.csv'))

for (fn in 1:NFACTORS){
  w_ordered<-w[, fn][order(abs(w[, fn]), decreasing = TRUE)]
  write.csv(w_ordered, paste0(nmf_outdir, '/top_weights/weights', fn, '.csv'))
  
  
}







#### Cluser ####
sel_factors=c(2,4)
sel_factors=c(1,2,3)

clusters_single <- kmeans(t(h)[,sel_factors], centers = 3)

covariates$cluster_s<-clusters_single$cluster[match(rownames(covariates),names(clusters_single$cluster))]
covariates$cluster_m<-clusters_mofa$cluster[match(rownames(covariates),names(clusters_mofa$cluster))]



chisq.test(clusters_single$cluster,covariates$COHORT )
chisq.test(clusters_single$cluster,covariates$COHORT )

df1=covariates
chisq.test(df1$cluster_m, df1$COHORT)

chisq.test(df1$cluster_s, df1$COHORT)

kruskal.test(df1$NHY, as.factor(df1$cluster_m ))
kruskal.test(df1$NHY, as.factor(df1$cluster_s ))

kruskal.test(df1$NP3_TOT, as.factor(df1$cluster_m ))
kruskal.test(df1$NP3_TOT, as.factor(df1$cluster_s ))

kruskal.test(df1$NP2_TOT, as.factor(df1$cluster_m ))
kruskal.test(df1$NP2_TOT, as.factor(df1$cluster_s ))




dim(df1)







