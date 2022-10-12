
library(ggplot2)
library(reshape2)
#install_phantomjs(version = "2.1.1",
#                  baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/",
#                  force = FALSE)
library('flextable')
#install.packages('webshot')
library('webshot')
library(dplyr)

library('readxl')
source('literature_review_combined.R')
sh1<-read_excel(paste0(data_int_dir,'MultiOmics Tools Evaluation_2.xlsx' ),sheet=1)
sh1<-sh1[,c(2,5:15)]
sh1<-sh1[-c(2),]


sh1 <- sh1[!((is.na(sh1[1]=='') | (sh1[1]=='') )),]
sh1<-data.frame(sh1)
rownames(sh1)<-sh1$Name


method_cats<-sh1

# remove empty rownames

colnames(method_cats)[1]<-'tool'
rownames(method_cats)<-method_cats$tool
method_cats$tool<-NULL
method_cats[method_cats=='X']<-1
method_cats[method_cats!=1]<-0
method_cats[is.na(method_cats)]<-0
# now melt 
method_cats_t<-  add_rownames(method_cats, var='rowname')


t2<-method_cats_t %>%
    gather(key, value, -rowname)%>%
    filter(value==1) 




t3<-t2 %>%  
  select(rowname, key)

t4<-plyr::ddply(t3, .(rowname), colwise(paste), collapse = ", ")

t4



dd<-method_cats
method_cats_ordered<-dd[
  with(dd, order(dd[,1], dd[,2], dd[,3], dd[,4], dd[,5],`Deep.Learning`,   decreasing = TRUE )),
]
method_cats_ordered
method_cats_ordered %>% as.data.frame() %>%
  add_rownames() %>% 
  flextable() 


sh2<-read_excel(paste0(data_int_dir,'MultiOmics Tools Evaluation_2.xlsx' ),sheet=2)

sh2<-sh2[c(2,5:10)]
sh2 <- sh2[!((is.na(sh2[1]=='') | (sh2[1]=='') )),]
sh2<-data.frame(sh2)
rownames(sh2)<-sh2$Name




sh2$Name<-NULL
sh2_chr<-sapply(sh2, as.character)
row.names(sh2_chr)<-row.names(sh2)
sh2$Ch1...Nonlinear.interactions[startsWith('X*',sh2$Ch1...Nonlinear.interactions)]
sh2_chr[which(startsWith(sh2_chr,'X', trim=TRUE))]<-1

sh2_chr[sh2_chr!=1]<-NA


sh2_chr %>% as.data.frame() %>%
  add_rownames() %>% 
  flextable() 

#EXPORT TO LATEX
sh2_chr_latex=sh2_chr
df<- sh2_chr_latex[,-c(5)]
df[is.na(df)]<-' '
colnames(df)<-c('Non-linear interactions', 'Unequal sizes', 'Missing Data', 'NP problem', 'Heterogeneous Datasets')
df[df==1]<-as.character('\\checkmark')
df<-as.data.frame(df)
df$x<-'\\\\'
write.table(df, file = "review/output/method_challenges.txt", sep=' & ', row.names = TRUE, quote = FALSE)

