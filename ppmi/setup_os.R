
os_dir='/Volumes/GoogleDrive/Other computers/My computer (1) (1)/'
if (Sys.info()['sysname']=='Darwin'){
  setwd(os_dir)
  data_dir<-os_dir
  
}else{
  setwd('D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/')
  data_dir<-'D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/'
  
  
}


