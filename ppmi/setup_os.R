
os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'
#setwd('/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/../')

isRStudio <- Sys.getenv("RSTUDIO") == "1"
if (isRStudio){
  script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
  script_dir<-paste0(script_dir, '/../')
  script_dir
}else{
  script_dir<- "D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/"
  
}


if (Sys.info()['sysname']=='Darwin'){
  setwd(os_dir)
  data_dir<-os_dir
  script_dir<-'/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/../'
  
}else if (Sys.info()['sysname']=='Linux'){
  setwd('/data8TB/efiath/git/mofa_tutorial/')
  data_dir<-'/data8TB/efiath/git/mofa_tutorial/'
  
}else  {
  setwd('D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/')
  data_dir<-'D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/'
  script_dir<-'D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/'
  
  
}


