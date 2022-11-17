

#install.packages('ggforce')
library(ggforce)
library("scales")
#install.packages('ggsci')
library("ggsci")



# TODO: change functionality, don't replace but add to more than one categories at the
rel_txt<-1.5



relabel_tools<-function(Var1){
  
  # This function groups instances of tools that are the same as others. 
  # The relabel has to exist in the matrix with tools and categories
  Var1<-sapply(Var1,function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*tcnn.*|.*classifier.*", 
                   ".*consensus*",
                   ".*autoencoder.*|.*pathcnn.*|.*umap.*|.*tsne.*|.*deep.*.|*neural.*|.*graph convolutional network.*|.*cdrscan.*|.*deepomix.*",
                   '.*pins.*|.*kernel.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*mofa.*',
                   '.*factor.*|.*decomposition.*|.*partitioning.*|.*pca.*|.*diverse.*', 
                   '.*netwo.*|.*piumet.*|.*omics integrator.*|.*inet.*|.*nem-tar.*',
                   '.*snf.*|.*coni.*|.*netdx.*|.*paradigm.*',
                   '.*nemo.*',
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*|.*panther.*|.*david.*|.*gsea.*|.*enrichment.*',
                   '.*metaboanalyst.*|.*adas.*|.*movics.*|.*mousse.*|.*timeg.*|.*miodin.*|.*ioda.*', 
                   '.*kmeans.*|.*k-means.*'), 
                 c( "ML/DL Classification", 
                    "ML Clustering", 
                    'pathcnn',
                    'cimlr',
                    'Regression',
                    'mofa',
                    'fa ', 
                    'interaction_nets',
                    'snf', 
                    'nemo',
                    'cca',
                    'Correlation',
                    'pls',
                    'multiomics pathway analysis',
                    'Other tools', 
                    'ML - DR'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(Var1)                 
}



group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*tcnn.*|.*classifier.*", 
                   ".*consensus*",
                   ".*autoencoder.*|.*pathcnn.*|.*umap.*|.*tsne.*|.*deep.*.|*neural.*|.*graph convolutional network.*|.*cdrscan.*|.*deepomix.*",
                   '.*pins.*|.*kernel.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*|.*partitioning.*|.*pca.*|.*diverse.*|.*pathme.*', 
                   '.*netwo.*|.*piumet.*|.*omics integrator.*|.*inet.*|.*nem-tar.*',
                   '.*snf.*|.*coni.*|.*netdx.*|.*paradigm.*',
                   '.*lra cluster*',
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*|.*panther.*|.*david.*|.*gsea.*|.*enrichment.*',
                   '.*metaboanalyst.*|.*nemo.*|.*adas.*|.*movics.*|.*mousse.*|.*timeg.*|.*miodin.*|.*ioda.*', 
                   '.*kmeans.*|.*k-means.*'), 
                 c( "ML/DL Classification", 
                    "ML Clustering", 
                    'JDR - NL',
                    'Kernel based',
                    'Regression', 
                    'Factor analysis ', 
                    'Network-Based',
                    'Network-Based - Similarity network', 
                    'Similarity based',
                    'JDR - LN - Correlation',
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
  #' 
  #'df[Var1]=as.factor(df[Var1])
  #' @param Var1: name of the column to replace
  #' @param df: name of the df 
  #' @return: the whole df 
  
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x), 
                 tolower(c( "ML/DL Classification", 
                    'JDR - NL',
                    'JDR - LN',
                    'Regression', 
                    'JDR - LN - Matrix Factorization', 
                    'Network-Based',
                    'Network-Based - Similarity network', 
                    'JDR - LN - Correlation',
                    'Correlation',
                    'JDR - LN - Partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 )),
                  c( "ML/DL Classification", 
                    'jDR - NL',
                    'jDR - LN',
                    'Regression', 
                    'jDR - MF', 
                    'NB',
                    'NB - SN', 
                    'jDR - CCA',
                    'Correlation',
                    'jDR - PLS',
                    'multiomics pathway analysis',
                    'Other tools'
                  ))}
  )
  #df[Var1]=as.factor(df[Var1])
  return(df)                 
}




group_methods_to_short_new<-function(df, Var1){
  #' Use only before plotting 
  #' 
  #'df[Var1]=as.factor(df[Var1])
  #' @param Var1: name of the column to replace
  #' @param df: name of the df 
  #' @return: the whole df 
  
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x), 
                 tolower(c( "ML/DL Classification", 
                            'JDR - NL',
                            'JDR - LN',
                            'Regression', 
                            'Matrix Factorization', 
                            'Network-Based',
                            'Network-Based - Similarity network', 
                            'Kernel-Based',
                            'Probability-based',
                            'JDR - LN - Correlation',
                            'Correlation',
                            'JDR - LN - Partial least squares',
                            'multiomics pathway analysis',
                            'Other tools'
                 )),
                 c( "ML/DL Classification", 
                    'jDR - NL',
                    'jDR - LN',
                    'Regression', 
                    'MF', 
                    'NB',
                    'NB - SN',
                    'KB',
                    'PR',
                    'jDR - CCA',
                    'COR',
                    'jDR - PLS',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #df[Var1]=as.factor(df[Var1])
  return(df)                 
}






relabel_objectives<-function(obj_col){
  
  obj_col <- mapvalues(obj_col, from=c("understand molecular mechanisms",
                                       "understand regulatory processes",   
                                 "extract patterns",
                                 "diagnosis/prognosis",
                                 "biomarker discovery",
                                 "si",
                                 "multiomics pathway analysis",
                                 
                                 "drug response prediction"), 
                     to=c("UN","RN", "CO", "DI", "BM", "SI","PA" ,"DR"))
    return(obj_col)
}



relabel_omics<-function(omics_col){
  
  omics_col <- mapvalues(omics_col, from=c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Metagenomics'), 
                       to=c('TR', 'GE', 'EP', 'PR', 'MB', 'MT'))
  return(obj_col)
}


group_objectives_method<-function(df, Var1){
  #'Group objective code column 
  #'These groups are for objective - method
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand molecular.*'),
                                  c('diagnosis/prognosis', 'understand molecular mechanisms')))
  return(df)
}


group_disease<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                c(".*alzheimer.*|.*amyotrophic.*|.*anxiety*|.*huntington.*|.*depressi.*|.*parkinson.*|.*autism.*|.*multiple sclerosis.*|.*epilepsy.*",
                  '.*cardio.*|.*heart.*|.*coronary.*|.*valve.*|.*atrial.*|.*stroke.*|.*blood pressure.*',
                  '.*bowel.*|.*hep.*|.*liver.*|.*nafld.*|.*crohn.*|.*ibd.*', 
                  '.*arthritis.*|.*osteo.*|.*facioscapulo.*', 
                  '.*diabetes.*', 
                  '.*pulmonar.*|.*lung.*|.*copd.*|.*smoking.*|.*pneumonia.*|.*cigarette.*|.*covid.*|.*traffic.*',
                  '.*metabolic.*|.*insulin.*|.*fasting.*',
                  '.*bladder.*|.*kidney.*|.*renal.*', 
                  '.*immune.*|.*sjögren.*|.*vogt.*|.*hiv.*|.*cytokine.*'),
                c('Nervous', 
                  'Cardiovascular',
                  'Gastrointestinal', 
                  'Musculoskeletal', 
                  'Endrocrine', 
                  'Pulmonary', 
                  'Metabolism', 
                  'Urinary', 
                  'Immune'))}
   
  )
  #new_col=as.factor(new_col)
  return(df)                 
}



### Only keep the most common combinations!! 
filter_common_groups<-function(df_by_group,top_n = c(8,10), x_group ){
  #' Returns the most common groups 
  #' freq_cutoff[1]: y-axis, freq_cutoff[2]: x-axis variables
  #' @return most_common_pairs, most_common_groups
  df_most_common<-df_by_group %>%
    group_by(Var1)  %>%
    summarise( Freq=sum(Freq) , .groups='drop') %>%
    as.data.frame()
  sel<-order(df_most_common$Freq, decreasing = TRUE)[1:top_n[1]]
  most_common_pairs<-df_most_common[,'Var1'][sel]
  
  df_most_common<-df_by_group %>%
    group_by_at(x_group)  %>%
    summarise( Freq=sum(Freq) , .groups='drop') %>%
    as.data.frame()
  sel<-order(df_most_common$Freq, decreasing = TRUE)[1:top_n[2]]
  most_common_groups<-df_most_common[,x_group][sel]
  
  return(list(most_common_pairs, most_common_groups))
}

 



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
#df<-df_to_plot

#[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
#[8] "#FF7F00"  "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
mycolors=c("epigenomics - transcriptomics" = '#3C5488FF', 
           "epigenomics - genomics" = '#4DBBD5FF',
           "genomics - transcriptomics" = '#91D1C2FF',
           "epigenomics - proteomics" = '#8491B4FF',
           "proteomics - transcriptomics" = '#00A087FF',
           "metabolomics - proteomics" = "#F39B7FFF" ,
           "metabolomics - transcriptomics" = "#B09C85FF", 
           "metabolomics - metagenomics" = "#7E6148FF", 
           "other" = "#808080")

mycolors=c("E - T" = '#3C5488FF', 
           "E - G" = '#4DBBD5FF',
           "G - T" = '#91D1C2FF',
           "E - P" = '#8491B4FF',
           "P - T" = '#00A087FF',
           "M - P" = "#F39B7FFF" ,
           "M - T" = "#B09C85FF", 
           "M - Mg" = "#7E6148FF", 
           "other" = "#808080")



plotbyObjective<-function(df, legend_t="Omics combinations", plot_width=8, 
                          plot_height=8, xangle=0, plot_cols=FALSE, stack_horiz=FALSE){ 
  text_size=12
  text_size_small=12
  fname=paste0('plots/barplot_byGroup', as.character(x_group), '_', colname,  
               '.jpeg')
  fname2=paste0('plots/barplot_byGroup', as.character(x_group), '_', colname,  
               '.jpeg')
  
  if (!('key_names' %in% colnames(df))){
    df['key_names']=df[x_group]
  }
  
  #[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
  #[8] "#FF7F00"  "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
  
  #### plot filter 
 
  df$y = df$Freq
  df$y = df$percent
  #g<-ggplot(df, aes(x=reorder(key_names, -y, sum), y=y, fill=Var1))+
  g<-ggplot(df, aes(x=reorder(key_names, -y, sum), y=y, fill=Var1))+
  
    geom_bar(stat='identity',position='stack', color='black')+
    guides(fill = guide_legend(title = legend_t))+
    labs(x=NULL, y='percentage of total omics combinations')
  
  if (stack_horiz){  g<-g+   coord_flip()}

  
  
  
  # text size related options
  #,
  g=g+ theme(axis.text.x = element_text(size=text_size, margin = margin(t = 0,b=0), vjust = 0.5,  angle = xangle,hjust=1))+
    theme(axis.text.y = element_text(size=text_size), axis.title.y =element_text(size=text_size),
    legend.text=element_text(size=text_size_small))
    
  
  if (colname == 'method'){
    print('paired brewer')
    g=g+ scale_fill_brewer(palette = 'Paired')+
      theme(plot.margin=unit(c(1,1,1,1.8),"cm"))
    
    #}else if (x_group=='disease_group'){
    #g=g+scale_fill_manual(values = mycolors)
    #g = g +theme(plot.margin=unit(c(1,1.8,4,2.2),"cm"))
  }else{
      #theme(legend.position = 'none')
    g=g+scale_fill_manual(values = mycolors)
      g = g +theme(plot.margin=unit(c(1,1,1,1.8),"cm"))
      

  }

    if ('Cancer' %in% colnames(df_to_plot)){

    if (plot_cols){
      print('add columns')
      # Add the two plots for cancer and not cancer in a row
      
      g<-g+facet_wrap(vars(Cancer),  scales = 'free_x',  labeller = labeller(Cancer=
             c('no'='Other Diseases','yes' ='Cancer')))+
        theme(strip.text.y = element_text(size = text_size))
      
        #scale_x_discrete(expand = c(0, 0.5))
      # g<-g+ facet_row(vars(Cancer), scales = 'free',space='free',labeller=
      #               labeller(Cancer=c('no'=' ','yes' =' '))) +
        
        
        # add the percentage labels 
      
    }else{
      print('add rows')
      
        g<-g+ facet_wrap( ~Cancer, ncol = 1,  labeller = labeller(Cancer=
                  c('no'='Other Diseases','yes' ='Cancer')),  
                       scales = 'free')   +
          theme(strip.text.x = element_text(size = text_size))
        
    }
      
      g<-g+ geom_text(aes( label=ifelse(percent>12,paste0(sprintf("%.0f", percent),"%"),'')), 
                      position=position_stack(vjust=0.5), size = 3) 
      
      ggsave(fname, width = plot_width, height=plot_height)
      ggsave(fname2, width = plot_width, height=plot_height)
      
      
    }
  

    return(g)
  
  
}


relabel_objectives_short<-function(df_to_plot){
  df_to_plot$labels<-df_to_plot$key_names
  ind<-df_to_plot$labels%in% c('extract patterns')
  df_to_plot[ind,]$labels<-'detect molecular patterns'
  df_to_plot[ind,]$key_names<-'detect molecular patterns'
  
  ind<-df_to_plot$labels%in% c('understand molecular mechanisms')
  df_to_plot[ind,]$labels<-'understand  molecular\n mechanisms'
  df_to_plot[ind,]$key_names<-'understand  molecular\n mechanisms'
  
  ind<-df_to_plot$labels%in% c('si')
  df_to_plot[ind,]$labels<-'subtype identification'
  df_to_plot[ind,]$key_names<-'subtype identification'
  
  ind<-df_to_plot$labels%in% c('understand regulatory processes')
  df_to_plot[ind,]$labels<-'understand regulatory\n processes'
  df_to_plot[ind,]$key_names<-'understand regulatory\n processes'
  
  return(df_to_plot)
}

relabel_omics_short<-function(x_data){
  #' Use only before plotting 
  #' 
  #'df[Var1]=as.factor(df[Var1])
  #' @param Var1: name of the column to replace
  #' @param df: name of the df 
  #' @return: the whole df 
  
  x_data<-sapply(x_data,function(x){
    mgsub::mgsub(tolower(x), 
                # metagenomics should be relabelled first
                   c("metabolomics - metagenomics",
                     "epigenomics - transcriptomics", 
                   "epigenomics - genomics",
                   "genomics - transcriptomics", 
                   "genomics - proteomics",
                   "epigenomics - proteomics" ,
                   "proteomics - transcriptomics", 
                   "metabolomics - proteomics" ,
                   "metabolomics - transcriptomics"),
                   c('M - Mg' ,
                  'E - T',
                    'E - G',
                    'G - T',
                    'G - P',
                    'E - P',
                    'P - T',
                    'M - P',
                    'M - T'
                    )
  )})
  return(x_data)                 
}




run_sankey<-function(df_to_plot,axis1, axis2,cancer_filter){
  

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
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_viridis_d()+
    ggtitle(paste0("Multi omics objectives, Cancer = ", cancer_filter))
  
  
  ggsave(paste0('plots/ggalluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png'), width = 7, height=6)
  
}



