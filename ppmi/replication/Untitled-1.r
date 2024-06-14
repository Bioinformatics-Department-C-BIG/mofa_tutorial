
data_dir
external_csv<-paste0(data_dir, '/ppmi/replication/csf-proteomics-mmc3_PROTEINS.csv')
external_data<-read.csv2( external_csv, sep = ',')

colnames(external_data)

markers<-c('VGF', 'SERPINF2', 'CHGA', 'SEPRINC1', 'PLG', 'FGA', 'FGB', 'APOC3', 'CFB', 'KNG1', 'GC', 'F2')

df_proteomics<-MOFAobjectPD@data[[3]]$group1
get_symbol_from_uniprot(rownames(df_proteomics))
rownames(df_proteomics)<-get_symbol_from_uniprot(rownames(df_proteomics))$SYMBOL

sm_pd<-samples_metadata(MOFAobjectPD)

to_Test<-df_proteomics[rownames(df_proteomics) %in% markers,]

cor_results<-corr.test(t(to_Test), sm_pd$NP3TOT_LOG)
cor_results$p.adj

intersect(external_data$PG.Genes, markers)



