


most_frequent_5
omic<-edge_list[1:2,]
top_combs=list()
top_combs_content=list()
for (i in 1:length(edge_list)){
  comb1<- stats %>% 
    filter(Cancer == cancer_filter) %>% 
    filter(Data %like%omic[i,]$X1  & Data %like% omic[i,]$X2  ) 
  print(omic[i,])
  print(comb1$Disease[!is.na(comb1$Disease)])
  print(comb1[,c('Title', 'ObjeMeth', 'Data', 'Citations')])
  top_combs_content[[i]]<-comb1[,c('Title', 'ObjeMeth', 'Data', 'Citations')]
  top_combs[[i]]<-comb1$Disease[!is.na(comb1$Disease)]
  top_combs
}
top_combs %>% count()
lapply(top_combs, write, "combinations_diseases.txt", append=TRUE, ncolumns=1000)

