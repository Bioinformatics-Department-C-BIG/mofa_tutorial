






s_stats<-read.csv( paste0(outdir_orig,'nmf_all_stats_mi.csv'))
colnames(s_stats)<-c('1','g_params', 'NFACTORS', 'k_centers', 'pvalue', 'cor', 'mut_inf', 'nf')
s_stats$g_params=as.factor(s_stats$g_params)

ggplot(s_stats, aes(x=k_centers, y=mut_inf))+
  geom_point(aes(x=k_centers, y=mut_inf, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~g_params)

  
