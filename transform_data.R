

highly_variable_genes_mofa[,1]
factor_vals =  get_factors(MOFAobject) 
factor_vals

MOFAobject@intercepts$mRNA

W_f1 = get_weights(MOFAobject)$mRNA[,'Factor1']
S1 = get_data(MOFAobject)$mRNA$group1[,1]
dotp = W_f1 %*%  S1  # Value for sample 1 
W_f1
S1
get_weights(MOFAobject)$mRNA$factor1 %*% get_data(MOFAobject)

get_expectations(MOFAobject)     
