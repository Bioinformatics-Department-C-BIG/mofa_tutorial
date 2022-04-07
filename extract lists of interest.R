


most_frequent_5
omic<-edge_list
top_combs=list()
for (i in 1:7){
  comb1<- stats %>% 
    filter(Cancer == cancer_filter) %>% 
    filter(Data %like%omic[i,]$X1  & Data %like% omic[i,]$X2  ) 
  print(omic[i,])
  print(comb1$Disease[!is.na(comb1$Disease)])
  top_combs[[i]]<-comb1$Disease[!is.na(comb1$Disease)]
  top_combs
}
top_combs %>% count()
lapply(top_combs, write, "combinations_diseases.txt", append=TRUE, ncolumns=1000)

