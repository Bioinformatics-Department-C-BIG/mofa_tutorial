
#################################### VISITS ##########################

#### load both visits to compare: 



merged_path_file

v08_paths<-read.csv(paste0( out_compare, 'V081p.adjust_FALSE', '.csv'))
bl_paths<-read.csv(paste0( out_compare, 'BL1p.adjust_FALSE', '.csv'))
dim(v08_paths); dim(bl_paths)

T_P<- 0.05
v08_paths_sig<-v08_paths[v08_paths$fish<T_P,]
bl_paths_sig<-bl_paths[bl_paths$fish<T_P,]
dim(v08_paths_sig)
dim(bl_paths_sig)
inter<-intersect(v08_paths_sig$Description,bl_paths_sig$Description )
#View(inter)  

inter

listInput_all_mods<-list(bl=bl_paths_sig$Description,
                         v08=v08_paths_sig$Description)



listInput<-listInput_all_mods
res_overlap<-calculate.overlap(listInput)
bl_only<-bl_paths_sig[!bl_paths_sig$Description %in% inter,]
v08_only<-v08_paths_sig[!v08_paths_sig$Description %in% inter,]
dim(v08_only)[1]


#View(bl_paths[bl_paths$Description %in% v08_only$Description,])
write.csv(bl_paths[bl_paths$Description %in% v08_only$Description,], 
          paste0(out_compare,'V08_only_single', T_P, '.csv'))




myCol <- brewer.pal(3, "Pastel2")[1:2]

venn.diagram(listInput,
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             cex=2.5,
             
             
             
             filename = paste0(out_compare,'all_modalities_','visits_venn_diagramm.png'), output=TRUE)




#### 



###get parents: 
library(GOfuncR)


### TODO: get the ids from my rna/prot files? 
### And get the parents or find where they are in the hierarchy
### is mofa giving me only upper? 
rna_parents<-get_parent_nodes(enrich_rna_sig$ID)
prot_parents<-get_parent_nodes(enrich_proteins_sig$ID)
mir_parents<-get_parent_nodes(enrich_mirnas_sig$ID)

rna_par_01<-rna_parents[rna_parents$distance==c(1),]
prot_parents_01<-prot_parents[prot_parents$distance==c(1),]
mir_parents_01<-prot_parents[mir_parents$distance==c(1),]

intersect(rna_par_01$child_go_id, prot_parents_01$child_go_id)
inter<-intersect(rna_par_01$parent_go_id, prot_parents_01$parent_go_id)
prot_parents_01[prot_parents_01$parent_go_id %in% inter,]

prot_parents_01
### we need the id, not the description 
rna_parents<-get_parent_nodes(merged_factors_mofa$Description)



### 

