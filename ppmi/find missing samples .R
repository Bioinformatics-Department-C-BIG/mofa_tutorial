

### find out why some mirnas measurements are not in metadata 

ind<-sapply(unique_s, 
            function(x){
              x2<-unlist(strsplit(x, '_'))[1]
              print(x2)
                grep(x2,colnames(mirnas_rpmmm))
            })
                    

rna_meta<-read.csv('ppmi/ppmi_data/rnaseq/meta_data.11192021.csv')
  
  
ind2<-unlist(ind)
all_without_meta<-colnames(mirnas_rpmmm)[ind2]

names_split<- strsplit(all_without_meta,split='\\.')
names_split_df<-do.call(rbind, names_split)

dim(names_split_df);
all_PATNO_wo_meta<-names_split_df[,3]
PATNO

MISSING<-rna_meta[gsub('PP-', '',rna_meta$PATNO) %in%all_PATNO_wo_meta,]

all_PATNO_wo_meta

### Are these MISSING samples in proteomics? 
## None of these samples are in proteomics so leave them for now unless we are doing rnaseq!! 

all_PATNO_wo_meta %in% ppmi_prot$PATNO
