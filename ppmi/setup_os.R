
os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'
#setwd('/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/../')

# isRStudio <- Sys.getenv("RSTUDIO") == "1"
# if (isRStudio){
#   script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#   script_dir<-paste0(script_dir, '/../')
#   script_dir
# }else{
#   script_dir<- "D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/"
  
# }


if (Sys.info()['sysname']=='Darwin'){
  setwd(os_dir)
  os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/mofa_tutorial/'
  os_dir='/Users/efiathieniti/Downloads/'

  script_dir<-'/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/../'
  ppmi_data<-'/Users/efiathieniti/Documents/ppmi_data/'
  ppmi_data<-'/Users/efiathieniti/Downloads/ppmi_data/'

  data_dir<-os_dir

  setwd(script_dir)
  
}else if (Sys.info()['sysname']=='Linux'){
  setwd('/data8TB/efiath/git/mofa_tutorial/')
  data_dir<-'/data8TB/efiath/git/mofa_tutorial/'
  script_dir<-'/data8TB/efiath/git/mofa_tutorial/'
  ppmi_data<-'/data8TB/efiath/git/mofa_tutorial/ppmi/ppmi_data/'
  
}else  {
  setwd('D:/DATADRIVE/Efi Athieniti/Documents/GitHub/mofa_tutorial/')
  data_dir<-'D:/DATADRIVE/Efi Athieniti/Documents/GitHub/mofa_tutorial/'
  ppmi_data<-'D:/DATADRIVE/Efi Athieniti/Documents/ppmi_backup/ppmi_data/'

  script_dir<-'D:/DATADRIVE/Efi Athieniti/Documents/GitHub/mofa_tutorial/'
  
  
}


