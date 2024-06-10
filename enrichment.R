

#install.packages("msigdbr")
library(msigdbr)
packageVersion("msigdbr")
packageVersion("clusterProfiler")

library(data.table)
#library(MOFAdata)
data("reactomeGS")

# GO annotations as present in the GO-basic obo file released on 2021-12-15 and NCBI gene2go annotations downloaded on 2022-01-03.


# Choose pathway gene list GO/KEGG etc 
#all_gene_sets = msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:BP')
category='C2';subcategory<-'CP:KEGG';
category='C5';subcategory<-'GO:MF';
category='C5';subcategory<-'GO:BP'; # total 7658 pathways



all_gene_sets = msigdbr(species = "Homo sapiens", category = category, subcategory =subcategory )
unique(all_gene_sets$gs_name) # total
df<-all_gene_sets[c('gs_name', 'ensembl_gene')]

df_wide<-data.table::dcast(as.data.table(df), 
                           gs_name ~ ensembl_gene)
df_wide<-as.data.frame(df_wide)
rownames(df_wide)<-df_wide$gs_name
df_wide$gs_name<-NULL


# make binary 
na_ind<-is.na(df_wide)
df_wide[na_ind]<-0
df_wide[!na_ind]<-1

df_wide<-as.matrix(df_wide)
df_wide %in% c(0,1)


#to_drop<-colnames(df_wide)[str_starts(colnames(df_wide), 'ASM')]
df_wide2<-df_wide[,-c(1:5)]


#### Some columns do not match... why? 

ind_to_keep<-!is.na(match(colnames(df_wide2), colnames(reactomeGS)))
length(ind_to_keep)
length(which(ind_to_keep))
df_wide3<-df_wide2[,ind_to_keep]
dim(df_wide3)

df_wide4<-as.matrix(apply(as.data.frame(df_wide3), 2,as.numeric))
rownames(df_wide4)<-rownames(df_wide3)

gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
write.csv(df_wide4,gs_file)




gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
gs_file_proteins<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')

gs_file_matrix<-read.csv(gs_file)
rownames(gs_file_matrix)

gs_file_matrix[,1]
colnames(gs_file_matrix)<-convert_to_gene_symbol(colnames(gs_file_matrix), view='RNA') 
write.csv(gs_file_matrix,gs_file_proteins, row.names=FALSE)


gs_file_matrix$X


