
out_compare<-'ppmi/plots/single/compare/'
enrich_compare<-paste0(out_compare, '/enrichment/')
dir.create(enrich_compare)

## TODO: compare proteins and mirnas 

####### Compare two time points ##### 
VISIT='V08'; process_mirnas=FALSE
source(paste0(script_dir, 'ppmi/config.R'))
gene_list_V08<-get_genelist_byVisit(VISIT)
VISIT='BL'; process_mirnas=FALSE; 
source(paste0(script_dir, 'ppmi/config.R'));outdir_s
gene_list_BL<-get_genelist_byVisit(VISIT)
VISIT='V04'
source(paste0(script_dir, 'ppmi/config.R'))
gene_list_V04<-get_genelist_byVisit(VISIT)


length(gene_list_BL)
length(gene_list_V08)

gse_compare<-compareCluster(geneClusters = list(BL=gene_list_BL,V08=gene_list_V08 ), 
              fun = "gseGO", 
               OrgDb='org.Hs.eg.db', 
               ont=ONT, 
               keyType = 'ENSEMBL') 

gse_compare
dot_comp<-dotplot(gse_compare, showCategory=30, split=".sign") + facet_grid(.~.sign)
dot_comp
ggsave(paste0(enrich_compare,prefix, 'dot_compare.jpeg' ), plot=dot_comp,
       dpi=300
)


gse_compare_x <- enrichplot::pairwise_termsim(gse_compare)
N_EMAP=200
emap_comp<-emapplot(gse_compare_x, showCategory=50,
                    cex.params = list(category_label = 1.1) ) 
emap_comp
ggsave(paste0(enrich_compare,prefix, 'emap_compare.jpeg' ), plot=emap_comp,
       dpi=300
       )


cnetplot(gse_compare)

#### to



