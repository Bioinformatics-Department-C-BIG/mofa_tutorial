

merged_targets

### Filter also by Differentially expressed genes only??? 
deseq2ResDF = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
outdir_s_mirnas<-outdir_s

process_mirnas=FALSE; source(paste0(script_dir, '/config.R'))
de_rnas = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)


merged_targets$symbol %in% outdir_s_mirnas
