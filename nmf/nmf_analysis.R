





#### Correlations ####


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
  png(paste0( nmf_param_str, '.png' ))
  
  corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")
  dev.off()  
  
  write.csv(stat[sig,], paste0('nmf/plots/cor_',mod, '.csv'  ))
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
  pheatmap::pheatmap(t(stat_filt), main="log10 adjusted p-values", cluster_rows = TRUE, color=col)
  
  
  
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
