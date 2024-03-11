proteins_outfile_plasma_raw = paste0(output_files, p_params_plasma , '.csv')


prot_untargeted_plasma_vsn_f
data<-read.csv2('/data8TB/efiath/git/mofa_tutorial//ppmi/output/proteomics_177_Plasma.csv', row.names=1)
head(data)

data_csf<-read.csv2('/data8TB/efiath/git/mofa_tutorial//ppmi/output/proteomics_177_Cerebrospinal Fluid.csv',row.names=1, check.names=FALSE)
data_csf['P00450',]
data['P00450',]
