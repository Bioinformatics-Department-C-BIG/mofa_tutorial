




# TODO: check that nrun==10!! and plot only the ones with 10 replicas
######## single NMF ####
mod='RNA'; selected_params<-'0.3_100_'; sel_NFACTORS=10
mod='miRNA'; selected_params<-'0.75_10_'; sel_NFACTORS=5


s_stats<-read.csv( paste0(outdir_orig,'nmf_all_stats_mi_', mod, '.csv'))
colnames(s_stats)<-c('1','nmf_params', 'NFACTORS', 'k_centers', 'pvalue', 'cor', 'mi', 'nf')
s_stats$nmf_params=as.factor(s_stats$nmf_params)
s_stats

### Plot mutual information 
ggplot(s_stats, aes(x=k_centers, y=mi))+
  geom_point(aes(x=k_centers, y=mi, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~nmf_params)+
  ylab('MI')


## Plot by nfactors 
ggplot(s_stats[s_stats$nmf_params==selected_params,], aes(x=NFACTORS, y=-log10(pvalue)))+
  geom_point(aes(x=NFACTORS, y=-log10(pvalue), color=k_centers))+ 
  geom_smooth()+
  facet_grid(~nmf_params)+
  ylab('-log10pvalue')


  
ggplot(s_stats[s_stats$nmf_params==selected_params,], aes(x=k_centers, y=-log10(pvalue)))+
  geom_point(aes(x=k_centers, y=-log10(pvalue), color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~nmf_params)+
  ylab('-log10pvalue')


# choose factors too 
ggplot(s_stats[s_stats$nmf_params==selected_params & s_stats$NFACTORS==sel_NFACTORS,], aes(x=k_centers, y=mi))+
  geom_point(aes(x=k_centers, y=mi, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~nmf_params)+
  ylab('MI')

# choose factors too 
ggplot(s_stats[s_stats$nmf_params==selected_params & s_stats$NFACTORS==sel_NFACTORS,], aes(x=k_centers, y=mi))+
  geom_point(aes(x=k_centers, y=mi, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~nmf_params)+
  ylab('MI')



########## MOFA ########


s_stats<-read.csv( paste0(outdir_orig,'all_stats_clusters.csv'), header = FALSE)
colnames(s_stats)<-c('1','TOP_PN', 'TOP_GN', 'MIN_COUNT_G', 'TOP_MN', 'MIN_COUNT_M', 'mofa_params', 'sel_coh_s','VISIT_S',  'scale_views',  'use_signif',
                            'run_mofa_complete', 'NFACTORS','cors_t' , 'max_cor','k_centers','mi' )


s_stats<-read.csv( paste0(outdir_orig,'all_stats_clusters.csv'))
colnames(s_stats)<-c('1','g_params', 'NFACTORS', 'k_centers', 'pvalue', 'cor', 'mut_inf', 'nf')
colnames(s_stats)<- c('1', 'TOP_PN', 'TOP_GN', 'MIN_COUNT_G', 'TOP_MN', 'MIN_COUNT_M', 'mofa_params', 'sel_coh_s','VISIT_S',  'scale_views',  'use_signif',
   'run_mofa_complete', 'N_FACTORS','cors_t' , 'max_cor','k_centers_m', 'mut_inf' )


#s_stats$g_params=as.factor(s_stats$g_params)

ggplot(s_stats, aes(x=k_centers_m, y=mi))+
  geom_point(aes(x=k_centers_m, y=mut_inf, color=N_FACTORS))+ 
  geom_smooth()+ 
  facet_grid(~TOP_GN)

ggplot(s_stats, aes(x=k_centers, y=mi))+
  geom_point(aes(x=k_centers, y=mi, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~TOP_GN)+
  ylab('MI')



ggplot(s_stats[s_stats$TOP_GN=='0.3' & !is.na(s_stats$NFACTORS) ,], aes(x=NFACTORS, y=mi))+
  geom_point(aes(x=NFACTORS, y=mi, color=k_centers))+ 
  geom_smooth()+
  facet_grid(~TOP_GN)+
  ylab('-log10pvalue')


# choose factors too 
ggplot(s_stats[s_stats$TOP_GN=='0.3'  & s_stats$NFACTORS==12,], aes(x=k_centers, y=mi))+
  geom_point(aes(x=k_centers, y=mi, color=NFACTORS))+ 
  geom_smooth()+
  facet_grid(~TOP_GN)+
  ylab('MI')




