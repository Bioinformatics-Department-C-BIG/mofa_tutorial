

# samples on columns 
metabolites<-read.csv('thal/ThalCypris_peak_metaboana,yst_2022_corrected_all_sample_cells - Copy.csv', row.names = 1)
metabolites
mir_data<-read.csv('thal/mature_counts.csv', row.names = 1)
mir_data<-t(mir_data)

### TODO: RUN VSN ONLY make a standalone function 
colnames(metabolites)<-gsub(".*\\.","",colnames(metabolites))

data_full<-list(miRNA=as.matrix(mir_data), 
                metabolites =as.matrix(metabolites))
colnames(metabolites)


# match 
match
MOFAobject <- create_mofa(data_full)

create_mofa(data_full)
