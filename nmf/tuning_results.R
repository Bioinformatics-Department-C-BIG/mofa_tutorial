




# TODO: check that nrun==10!! and plot only the ones with 10 replicas

s_stats<-read.csv( paste0(outdir_orig,'nmf_all_stats_mi.csv'))
colnames(s_stats)<-c('1','g_params', 'NFACTORS', 'k_centers', 'pvalue', 'cor', 'mut_inf', 'nf')
s_stats$g_params=as.factor(s_stats$g_params)


### Plot mutual information 
ggplot(s_stats, aes(x=k_centers, y=mut_inf))+
  geom_point(aes(x=k_centers, y=mut_inf, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~g_params)


### Plot mutual information 
ggplot(s_stats, aes(x=k_centers, y=cor))+
  geom_point(aes(x=k_centers, y=mut_inf, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~g_params) 









s_stats<-read.csv( paste0(outdir_orig,'all_stats_clusters.csv'))
colnames(s_stats)<-c('1','g_params', 'NFACTORS', 'k_centers', 'pvalue', 'cor', 'mut_inf', 'nf')
colnames(s_stats)<- c('1', 'TOP_PN', 'TOP_GN', 'MIN_COUNT_G', 'TOP_MN', 'MIN_COUNT_M', 'mofa_params', 'sel_coh_s','VISIT_S',  'scale_views',  'use_signif',
   'run_mofa_complete', 'N_FACTORS','cors_t' , 'max_cor','k_centers_m', 'mut_inf' )


#s_stats$g_params=as.factor(s_stats$g_params)

ggplot(s_stats, aes(x=k_centers_m, y=mut_inf))+
  geom_point(aes(x=k_centers_m, y=mut_inf, color=N_FACTORS))+ 
  geom_smooth()+ 
  facet_grid(~TOP_GN)




