



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


#### Compare also by group 


gse_compare<-compareCluster(geneClusters = list(BL=gene_list_BL,V08=gene_list_V08 ), 
              fun = "gseGO", 
               OrgDb='org.Hs.eg.db', 
               ont=ONT, 
               keyType = 'ENSEMBL') 

enrich_compare_path=paste0(enrich_compare,prefix)


plot_enrich_compare(gse_compare,enrich_compare_path)


#### to



