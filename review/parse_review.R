
install.packages("textreadr")
install.packages("tidyverse")
install.packages("pdfsearch")
library(textreadr)
library(tidyverse)
library(pdfsearch)
library(data.table)

dirct <- '/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Multi omics/automated/'
result <- keyword_directory(dirct, 
#keyword = 'proteom|transcriptom|RNA|miRNA|lipidomic|metabolomic|16|microbiome|metagenomic|genomic|DNA seq|genetic',
keyword = 'objective|aim',
surround_lines = 0, full_names = TRUE,max_search = 5)
head(result[,c('ID','line_text')], n = 100)
result$line_text2<-unlist(result$line_text)
result2<-result[,c('ID','pdf_name','line_text2')]
write.csv2(result2, 'results3.csv')
