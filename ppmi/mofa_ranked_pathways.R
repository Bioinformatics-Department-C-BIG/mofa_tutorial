
## ## TODO: if enrichment is already run then just load results
## load res. positive
# load cohort cors


#### choose factors from MOFA and RANK the pathways accordingly 
cohort_cors<-cors_pearson[,'CONCOHORT'] # TODO: LOAD or recalc
dim(res.positive$pval.adj)
length(cohort_cors)
var_captured  <-vars_by_factor_all$r2_total[[1]] # TODO: LOAD or recalc
vars_by_factor_all$r2_total
cor_t<-0.1
cor_t<-0.17

cohort_cors

sel_factors<-c('Factor3', 'Factor4', 'Factor5')
sel_factors<-c('Factor3', 'Factor4', 'Factor5')
sel_factors<-c( 'Factor3', 'Factor4')



sel_factors<-c('Factor3', 'Factor4', 'Factor7')

sel_factors<-which(abs(cohort_cors)>cor_t)


sel_factors<-which(abs(cohort_cors)>cor_t)
sel_factors
paste0(sel_factors)
cohort_cors[sel_factors]
var_captured
sel_factors<-which(abs(cohort_cors)>cor_t)
sel_factors
### TODO: add var captured
cohort_cors[sel_factors]

pos=TRUE
if (pos){
  enrich_weights<-res.positive$pval.adj # todo: load res.positive 
}else{
  enrich_weights<-res.negative$pval.adj
  
}
log10pvals<--log10(enrich_weights[, sel_factors])
abs(cohort_cors[sel_factors])
log10pvals[1]*abs(cohort_cors[sel_factors])[1]

log10pvals[1]
hist(log2(abs(cohort_cors[sel_factors])))

get_weighted_pvals<-function(enrich_weights){
  weighted_pvals<--log10(enrich_weights[, sel_factors])*abs(cohort_cors[sel_factors])
  head(weighted_pvals)
  pvals_mofa<-enrich_weights[, sel_factors]
  pvals_mofa_melted<-melt(pvals_mofa,value.name='p.adj')
  
  
  melted<-melt(weighted_pvals,value.name='weighted')
  head(pvals_mofa_melted)
  head(melted)
  melted_merged<-merge(pvals_mofa_melted,melted, by=c('Var1', 'Var2') )
  
  head(melted_merged)
  log2(0.2)
  ### RANKING
  melted_ord<-melted_merged[order(abs(melted_merged$weighted), decreasing=TRUE),]
  T=0.01
  melted_ord_sig<-melted_ord[melted_ord$p.adj<T,]
  
  dim(melted_ord_sig)
  length(unique(melted_ord_sig$Var1))
  
  write.csv(melted_ord_sig,paste0(out_compare, 'mofa_',pos,T, '.csv' ) )
  return(melted_ord_sig)
}

###   TODO: load res.positive from file



mofa_file<-paste0(out_compare, 'mofa_',T, cor_t , mofa_params, TISSUE) 
pos_ord<-get_weighted_pvals(res.positive$pval.adj)
neg_ord<-get_weighted_pvals(res.negative$pval.adj)

all_ord<-rbind(neg_ord,pos_ord )
all_ord<-all_ord[order(all_ord$weighted),]
### melted_ord<-melted_merged[order(-abs(melted_merged$p.adj), decreasing=TRUE),]

all_ord<-all_ord[order(all_ord$weighted, decreasing = TRUE),]


all_ord$Var1<-gsub('GOBP_', '', all_ord$Var1)
all_ord$Var1<-tolower(gsub('_', ' ', all_ord$Var1))
length(unique(all_ord$Var1))
colnames(all_ord)[1]<-'Description'
colnames(all_ord)

#mofa_enrich_file<-paste0(outdir,'/enrichment/', 'ranked_list', cor_t, '.csv')

mofa_enrich_dir=paste0(outdir,'/enrichment/', mofa_params, TISSUE, 'ranked_list', cor_t)

write.csv(all_ord, paste0(mofa_enrich_dir,'.csv'), row.names = FALSE)
Npaths<-25
all_ord$p.adj
all_ord$log10padj<--log10(all_ord$p.adj)

mofa_enrich_plot<-ggplot(all_ord[1:Npaths, ], 
                         aes( x=reorder(Description, weighted),
                                                         y=log10padj))+
geom_bar(position='dodge', stat='identity', width=0.7/4, fill='darkgreen')+
  theme(axis.title.y=element_blank(), 
        axis.text.y= element_text(size=15,color='black' ), 
        axis.text.x= element_text(size=15,color='black' ),
        legend.title =element_blank(), 
        legend.text = element_text(size=15))+
  labs(y='-log10pvalue')+
  scale_fill_gradient()+
 scale_fill_discrete(label=c('MOFA'))+

  coord_flip()
mofa_enrich_plot

##subtitle <- textGrob("subtitle", x=0, hjust=0, gp=gpar( fontface="italic"))
## TODO: add legends
ggsave(paste(mofa_enrich_dir, '.jpeg'),mofa_enrich_plot, dpi=300,
     width=Npaths/2+1,height=Npaths/2.5 )

ggsave(paste(mofa_file, '.jpeg'),mofa_enrich_plot, dpi=300,
       width=Npaths/2+1,height=Npaths/2.5 )

mofa_enrich_plot

unique(all_ord$Description)




