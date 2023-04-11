



if (use_signif){
  TOP_GN=0.50
  TOP_MN=0.75
}else{
  TOP_GN=0.10
  TOP_MN=0.5
  TOP_MN=0.75
}



### MOFA definitions
highly_variable_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa.csv')
highly_variable_sign_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa_signif.csv')

