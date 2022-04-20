




group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*tcnn.*|.*classifier.*",  
                   ".*autoencoder.*|.*pathcnn.*|.*umap.*|.*tsne.*|.*deep.*.|*neural.*|.*graph convolutional network.*|.*cdrscan.*|.*deepomix.*",
                   '.*cluster.*|.*pins.*|.*kernel.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*|.*partitioning.*|.*pca.*|.*diverse.*|.*pathme.*', 
                   '.*netwo.*|.*piumet.*|.*omics integrator.*|.*inet.*|.*nem-tar.*',
                   '.*snf.*|.*coni.*|.*netdx.*|.*paradigm.*',
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*|.*panther.*|.*david.*|.*gsea.*|.*enrichment.*',
                   '.*metaboanalyst.*|.*nemo.*|.*adas.*|.*movics.*|.*mousse.*|.*timeg.*|.*miodin.*|.*ioda.*', 
                   '.*kmeans.*|.*k-means.*'), 
                 c( "ML Classification", 
                    'JDR - NL',
                    'JDR - LN',
                    'Regression', 
                    'JDR - LN - Matrix Factorization', 
                    'Network-Based ',
                    'Network-Based - Similarity network', 
                    'JDR - correlation',
                    'Correlation',
                    'JDR - LN - Partial least squares',
                    'multiomics pathway analysis',
                    'Other tools', 
                    'ML - DR'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}

group_methods_to_short<-function(df, Var1){
  #' Use only before plotting 
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x), 
                 tolower(c( "ML Classification", 
                    'JDR - NL',
                    'JDR - LN',
                    'Regression', 
                    'JDR - LN - Matrix Factorization', 
                    'Network-Based ',
                    'Network-Based - Similarity network', 
                    'JDR - correlation',
                    'Correlation',
                    'JDR - LN - Partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 )),
                  c( "ML Classification", 
                    'JDR - NL',
                    'JDR - LN',
                    'Regression', 
                    'JDR - LN - MF', 
                    'NB ',
                    'NB - SNF', 
                    'JDR - correlation',
                    'Correlation',
                    'JDR - LN - PLS',
                    'multiomics pathway analysis',
                    'Other tools'
                  ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}







relabel_objectives<-function(obj_col){
  
  obj_col <- mapvalues(obj_col, from=c("understand molecular mechanisms",
                                       "understand regulatory processes",   
                                 "connect molecular patterns to phenotypic traits",
                                 "diagnosis/prognosis",
                                 "biomarker discovery",
                                 
                                 "subtype identification",
                                 "multiomics pathway analysis",
                                 
                                 "drug response prediction"), 
                     to=c("UN","UN", "CO", "DI", "BM", "SI","PA" ,  "DR"))
    return(obj_col)
}




group_objectives_method<-function(df, Var1){
  #'Group objective code column 
  #'These groups are for objective - method
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand.*'),
                                  c('diagnosis/prognosis', 'understand molecular mechanisms')))
  return(df)
}


group_disease<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                c(".*alzheimer.*|.*amyotrophic.*|.*anxiety*|.*depressi.*|.*parkinson.*|.*autism.*|.*multiple sclerosis.*|.*epilepsy.*",
                  '.*cardio.*|.*heart.*|.*coronary.*|.*valve.*|.*atrial.*',
                  '.*bowel.*|.*hep.*|.*liver.*|.*nafld.*|.*crohn.*', 
                  '.*arthritis.*|.*osteo.*|.*fasioscapulo.*', 
                  '.*diabetes.*', 
                  '.*pulmonar.*|.*lung.*|.*copd.*|.*smoking.*|.*cigarette.*|.*traffic.*',
                  '.*metabolic.*|.*insulin.*|.*fasting.*',
                  '.*bladder.*|.*kidney.*|.*renal.*'),
                c('Nervous', 
                  'Cardiovascular',
                  'Gastrointestinal', 
                  'Musculoskeletal', 
                  'Endrocrine', 
                  'Pulmonary', 
                  'Metabolism', 
                  'Urinary' ))}
   
  )
  #new_col=as.factor(new_col)
  return(df)                 
}



### Only keep the most common combinations!! 
filter_common_groups<-function(df_by_group,freq_cutoff = c(17,17) ){
  #' Returns the most common groups 
  #' freq_cutoff[1]: y-axis, freq_cutoff[2]: x-axis variables
  df_most_common<-df_by_group %>%
    group_by(Var1, Cancer)  %>%
    filter( sum(Freq) >= freq_cutoff[1]) %>%
    group_by_at(x_group)  %>%
    filter( sum(Freq) >= freq_cutoff[2])
  return(df_most_common)
}

 
library(ggforce)



adjust_facet_width<-function(g,fname){# convert ggplot object to grob object

  print('adjusting')
  gp <- ggplotGrob(g)
  
  # optional: take a look at the grob object's layout
  gtable::gtable_show_layout(gp)
  
  # get gtable columns corresponding to the facets (5 & 9, in this case)
  facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
  
  # get the number of unique x-axis values per facet (1 & 3, in this case)
  x.var <- sapply(ggplot_build(g)$layout$panel_scales_x,
                  function(l) length(l$range$range))
  
  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
  
  
  png(fname, height =plot_height*100, width=plot_width*100)
  # plot result
  grid::grid.draw(gp)
  dev.off()
  gg<-arrangeGrob(gp)
  ggsave('tmp_gg_save.png', height =plot_height, width=plot_width)
  
  return(gg)
}

plot_width=8
plot_height=8
angle=30
plot_cols=TRUE
df<-df_to_plot
plotbyObjective<-function(df, legend_t="Omics combinations", plot_width=8, plot_height=8, angle=30, plot_cols=FALSE){ 
  

  
  fname=paste0('plots/barplot_byGroup', as.character(x_group), '_', colname,  
               '.png')
  
  if (!('key_names' %in% df)){
    df['key_names']=df[x_group]
  }
  
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(15)

  #### plot filter 
  
  g<-ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack', color='black')+
    scale_fill_brewer(palette = 'Paired')+
    #scale_fill_manual(mycolors)+
    guides(fill = guide_legend(title = legend_t))+
    
    labs(x=NULL)+
    
    theme(axis.text.x = element_text(size=rel(1.5),angle = angle, vjust = 0.5, hjust=1))+
    theme(axis.text.y = element_text(size=rel(1.5)))+
    
    theme(plot.margin=unit(c(1,1,1.8,2.2),"cm"))+
    theme(legend.text=element_text(size=rel(1.4)))
   
    if ('Cancer' %in% colnames(df_to_plot)){

    if (plot_cols){
      print('add columns')
      
      #g<-g+facet_wrap(vars(Cancer),  scales = 'free',  labeller = labeller(Cancer=
             #c('no'=' ','yes' =' ')))+
       #scale_x_discrete(expand = c(0, 0.5))
      g<-g+ facet_row(vars(Cancer), scales = 'free',space='free',labeller=
                      labeller(Cancer=c('no'=' ','yes' =' ')))                                                                             

      ggsave(fname, width = plot_width, height=plot_height)
      
    }else{
      print('add rows')
      
        g<-g+ facet_wrap( ~Cancer, ncol = 1,  labeller = labeller(Cancer=
                  c('no'='Other Diseases','yes' ='Cancer')),  
                       scales = 'free_y')    }
      ggsave(fname, width = plot_width, height=plot_height)
      
      
  }
  
  


    return(g)
  
  
}


relabel_objectives_short<-function(df_to_plot){
  df_to_plot$labels<-df_to_plot$key_names
  ind<-df_to_plot$labels%in% c('connect molecular patterns to phenotypic traits')
  df_to_plot[ind,]$labels<-'connect molecular patterns to \n phenotypic traits'
  df_to_plot[ind,]$key_names<-'connect molecular patterns to \n phenotypic traits'
  
  ind<-df_to_plot$labels%in% c('understand molecular mechanisms')
  df_to_plot[ind,]$labels<-'understand \n molecular mechanisms'
  df_to_plot[ind,]$key_names<-'understand \n molecular mechanisms'
  
  
  return(df_to_plot)
}


run_sankey<-function(df_to_plot,axis1, axis2,cancer_filter, col){
  
  
  
  
  
  counts<-df_to_plot %>%
    ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))
  
  
  df<-counts
  
  ggplot(data = df, aes(x = type, stratum = stratum, alluvium = alluvium, y = n)) +
    # geom_lode(width = 1/6) +
    geom_flow(aes(fill = col), width = 1/6, color = "darkgray",
              curve_type = "cubic") +
    # geom_alluvium(aes(fill = stratum)) +
    geom_stratum(color = "grey", width = 1/6) + 
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    theme(
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_viridis_d()+
    ggtitle(paste0("Multi omics objectives, Cancer = ", cancer_filter))
  
  
  ggsave(paste0('plots/ggalluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png'), width = 7, height=6)
  
}

