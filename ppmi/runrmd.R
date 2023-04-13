
output_dir<-outdir_s
output_file<-paste0(outdir_s, '/deseq_results_static.html')

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))
print(script_dir)

rmarkdown::render(paste0(script_dir,'/deseq_results_static.Rmd'),
                  output_file = output_file,
                  output_dir =outdir_s, 
                     )



                  