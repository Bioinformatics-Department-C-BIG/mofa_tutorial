##### EMAP utils #### 


remove_subcomponents<-function(g, subcomp_min_edge=2){
  #'
  #'
  #'
  sub_gs<-components(g)$membership
  small_sub <- names(which(table(sub_gs) <= subcomp_min_edge))
  #get names of nodes to rm
  rm_nodes <- which(sub_gs %in% small_sub)
  #remove nodes by name
  g_filt<- delete_vertices(g, rm_nodes)
  return(list(g_filt, rm_nodes))
}
