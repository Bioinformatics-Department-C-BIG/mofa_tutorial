

formula_deseq_format='n'
cell_corr_deseq = TRUE
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



