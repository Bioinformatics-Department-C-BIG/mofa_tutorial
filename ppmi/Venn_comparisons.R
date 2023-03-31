

#### Metascripts 
#install.packages('UpSetR')
library('UpSetR')
library('dplyr')
#install.packages('VennDiagram')
library('VennDiagram')
library(grid)

### Table of samples from all visits 

out_compare<-'ppmi/plots/single/compare/'
dir.create(out_compare)


#### Firstly print numbers of samples 
visits<-c('BL', 'V04', 'V06', 'V08')
all_vs_ps<-filter_se(se,visits, sel_coh )

meta<-colData(all_vs_ps)
get_stats<-meta[,c('PATNO', 'EVENT_ID')]
table_counts<-table(get_stats)
common_patients<-rownames(table_counts[rowSums(table_counts)==dim(table_counts)[2],])


list_all_vs<-split(get_stats, f = get_stats$EVENT_ID)



ns<-lapply(list_all_vs, function(x){
  length(unique(x$PATNO))
})
patient_lists<-lapply(list_all_vs, function(x){
  x$PATNO
})
upset(fromList(patient_lists))
counts_table<-table(get_stats)

## TODO: which are the common samples in all the lists


#View(c(meta$PATNO,meta$EVENT_ID))
library(plyr)
meta %>% 
  filter(EVENT_ID==visits[1])



dirs<-dir(path = "ppmi/plots/single/", pattern = 'mirnas', 
    full.names = TRUE)
visits=c('BL', 'V04', 'V06', 'V08')

if  (process_mirnas){
  outdir_all<-paste0(outdir_orig, '/single/', 'mirnas_',visits, '_', m_params ,'coh_', sel_coh_s, '_',des)
  title='mirnas'
}else{
  outdir_all<-paste0(outdir_orig, '/single/',  'rnas_',visits, '_', g_params ,'coh_', sel_coh_s, '_',des)
  title='rnas'
  
}


all_visits<-lapply(outdir_all, function(x) {
  outfile<-paste0(x, '/significant.csv')
  as.data.frame(read.csv(outfile))
  })

names(all_visits)
all_visits
list_of_mirs<-lapply(all_visits,function(df)
  df[df$sign_lfc=='Significant',]$X
)




listInput <- list(BL = list_of_mirs[[1]],
                  V04 =  list_of_mirs[[2]], 
                  V06 =  list_of_mirs[[3]],
                  V08 =  list_of_mirs[[4]])


jpeg(paste0(out_compare, 'upSet_diagram.jpeg'), res=200, width=800, height=500)
up<-upset(fromList(listInput))
up+grid.text('mirnas')
dev.off()

#data_with_intersection <- listInput %>%
#  unite(col = "intersection", -c("entry"), sep = "")


getwd()
venn.diagram(listInput,   
             filename = paste0(out_compare,'14_venn_diagramm.png'), output=TRUE)


