
library(ggplot2)
library(reshape2)
#install_phantomjs(version = "2.1.1",
#                  baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/",
#                  force = FALSE)
library('flextable')
#install.packages('webshot')
library('webshot')
library(dplyr)
method_cats<-read.csv2('review/datasets/method_categories_raw.csv', sep=',')

# remove empty rownames
method_cats <- method_cats[!((is.na(method_cats[1]=='') | (method_cats[1]=='') )),]
colnames(method_cats)[1]<-'tool'
rownames(method_cats)<-method_cats$tool
method_cats$tool<-NULL
method_cats[method_cats=='X']<-1
method_cats[method_cats!=1]<-0



dd<-method_cats
method_cats_ordered<-dd[
  with(dd, order(dd[,1], dd[,2], dd[,3], dd[,4], dd[,5], dd[,6], dd[,7], decreasing = TRUE )),
]
method_cats_ordered
method_cats_ordered %>% as.data.frame() %>%
  add_rownames() %>% 
  flextable() 






