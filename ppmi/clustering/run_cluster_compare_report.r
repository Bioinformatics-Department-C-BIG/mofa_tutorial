

formula_deseq_format='n'
cell_corr_deseq = FALSE
process_mirnas = FALSE

for (VISIT_COMP in c('V08', 'V06', 'BL', 'V04')){

        rmarkdown::render("ppmi/clustering/cluster_compare_report.Rmd", params = list(
            title=paste0('Mofa Clusters comparisons - ', VISIT_COMP, ' - ', cell_corr_deseq, ' - ', formula_deseq_format ), 
            VISIT_COMP=VISIT_COMP,
            cell_corr_deseq=cell_corr_deseq,
            formula_deseq_format=formula_deseq_format

        ), 
        output_file = paste0("Cluster comparisons-", VISIT_COMP, "-", cell_corr_deseq, '-',formula_deseq_format, ".html")

        )
}




## All timepoints per cluster 
## all modalities - only Volcano plots because we dont have enrichment for proteins yet?

formula_deseq_format='all'
for (VISIT_COMP in c('V08', 'V06', 'BL', 'V04')){
#for (VISIT_COMP in c('V08', 'V06')){

    # 
    rmarkdown::render("ppmi/clustering/cluster_compare_report_all_modalities.Rmd", params = list(
            title=paste0('Mofa Clusters comparisons - all modalities ', VISIT_COMP, ' - ', cell_corr_deseq, ' - ', formula_deseq_format ), 
            VISIT_COMP=VISIT_COMP,
            cell_corr_deseq=cell_corr_deseq,
            formula_deseq_format=formula_deseq_format

        ), 
        output_file = paste0("Cluster comparisons-all modalities ", VISIT_COMP, "-", cell_corr_deseq, '-',formula_deseq_format, ".html")

        )
}



## Plot annotation of factors 
clust_metric<-'NP2PTOT_LOG'
clust_metric<-'updrs3_score_on'
nf<-get_factors_for_metric(clust_metric)
 #nf = c(1,11,12) # factors to get enrichment 
 nf_s<-paste0( nf, collapse = ' ')

 rmarkdown::render("ppmi/mofa_enrich_report.Rmd", params = list(
            title=paste0('MOFA enrich report - factors annotation ',clust_metric, ', ', nf_s  ), 
            nf = nf),
                    output_file = paste0("Mofa Enrichment, ", clust_metric, ', ', nf_s, ".html")
 )







for (VISIT_COMP in c('V08', 'V06', 'BL', 'V04')){

        rmarkdown::render("ppmi/clustering/cluster_compare_report.Rmd", params = list(
            title=paste0('Mofa Clusters comparisons - ', VISIT_COMP, ' - ', cell_corr_deseq, ' - ', formula_deseq_format ), 
            VISIT_COMP=VISIT_COMP,
            cell_corr_deseq=cell_corr_deseq,
            formula_deseq_format=formula_deseq_format

        ), 
        output_file = paste0("Cluster comparisons-", VISIT_COMP, "-", cell_corr_deseq, '-',formula_deseq_format, ".html")

        )
}

# RUN MOFA ENRICH REPORT
Sys.setenv(RSTUDIO_PANDOC="--- insert directory here ---")
Sys.getenv("RSTUDIO_PANDOC")

Sys.setenv(RSTUDIO_PANDOC="--- insert directory here ---")

        rmarkdown::render("ppmi/mofa_enrich_report.Rmd", params = list(
      
                outdir = outdir,
                title= 'Mofa ', 
                nf= c(13,15,20,24,28)

        ), 
        output_file = paste0("mofa enrich report.html")

        )
