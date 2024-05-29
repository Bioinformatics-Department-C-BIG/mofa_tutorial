



genes_files<-list.files(paste0(paste0(outdir, '/enrichment/pcgse/')), pattern = 'feature.*pos_neg.csv',full.names = TRUE)
genes_files


genes_all_modes= list()
all_Fs<-list()

mod1 = convert_to_gene_symbol(concatenate_top_features(MOFAobject,factors_all=c(24), view = 'RNA', top_fr=0.01)$feature, view='RNA')
mod2 =convert_to_gene_symbol(concatenate_top_features(MOFAobject,factors_all=c(24), view = 'proteomics_t_plasma', top_fr=0.05)$feature,  view ='proteomics_t_plasma' )
mod3 = convert_to_gene_symbol(concatenate_top_features(MOFAobject,factors_all=c(24), view = 'proteomics_t_csf', top_fr=0.05)$feature,  view ='proteomics_t_csf' )
mod4 = convert_to_gene_symbol(concatenate_top_features(MOFAobject,factors_all=c(24), view = 'proteomics_csf', top_fr=0.05)$feature, view = 'proteomics_csf')
write.csv(c(mod1, mod2, mod3, mod4),paste0(outdir, '/top_weights/top_weights_all_mods.csv'))


mod4_u = concatenate_top_features(MOFAobject,factors_all=c(24), view = 'proteomics_csf', top_fr=0.05)$feature

convert_to_gene_symbol(mod4_u, view = 'proteomics_csf' , conv_uniprot=TRUE)

for (genes_file in genes_files ) {

    name_list<-unlist(regmatches(genes_file,regexec('features_.*_GO',genes_file )))

    genes_all_modes[[name_list]]<-as.data.frame(read.csv(genes_file))
    genes_all_modes[[name_list]]$name=name_list
    print(head(genes_all_modes[[name_list]]))


    



}
all_fs<-do.call(rbind,genes_all_modes)



factor = 24

for (factor in 1:34){

        only_1f<-all_fs[all_fs$ind == factor,]
        only_1f$values



        #BiocManager::install('STRINGdb')
        library(STRINGdb)
        #diff_exp_example1 = read.csv('/tmp/Rtmpm1GKNE/vscode-R/F11')
        diff_exp_example1 = only_1f

        string_db <- STRINGdb$new( version="12.0", species=9606, score_threshold=200, network_type="full", input_directory="")
        example1_mapped <- string_db$map( diff_exp_example1, "values", removeUnmappedRows = TRUE )


        hits <- example1_mapped$STRING_id

        png(paste0(outdir,'/string', factor, '.png'), res=300, units='in',  width=10, height=10)
        string_db$plot_network( hits , required_score = 0.9)


        graphics.off()      
}

# get enrichment analysis 
enrichment <- string_db$get_enrichment( hits )
colnames(enrichment)

unique(enrichment$category)
enrichment_new<-enrichment[enrichment$p_value<0.000000001,]
enrichment_new<-enrichment_new %>% filter(category!='PMID')
head(enrichment_new[,c('category', 'description')], 200)
#enrichment_new$inputGenes


example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, grepl('proteomics_csf',diff_exp_example1$name)),
 logFcColStr='ind' )


    payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color )

png(paste0(outdir,'/string', factor, '.png'), res=300, units='in',  width=10, height=10)

    string_db$plot_network( hits , required_score = 0.7,  payload_id=payload_id)
dev.off()








