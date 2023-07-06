
plot_top_weights2 <- function(object, view = 1, factors = 1,
                             nfeatures = 10, abs = TRUE, scale = TRUE, sign = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (nfeatures <= 0) stop("'nfeatures' has to be greater than 0")
  if (sign=="all") { abs <- TRUE}
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  
  # Get views
 # view <- .check_and_get_views(object, view)
  
  # Get factor names
  #factors <- .check_and_get_factors(object, factors)
  #scale=TRUE
  # Collect expectations  
 # W<-get_weights(MOFAobject_gs,factors=1, views = 'RNA', scale = TRUE, as.data.frame=TRUE)
  
  W <- get_weights(object, factors = factors, views = view, as.data.frame=TRUE)
  
  
 # W$value
  
  # Scale values by weight with highest (absolute) value
  if (scale) W$value <- W$value/max(abs(W$value))
  
  # Store sign
  W <- W[W$value!=0,]
  W$sign <- ifelse(W$value>0, "+", "-")
  
  # Select subset of only positive or negative weights
  #if (sign=="positive") { W <- W[W$value>0,] } else if (sign=="negative") { W <- W[W$value<0,] }
  
  # Absolute value
  if (abs) W$value <- abs(W$value)
  
  # Extract relevant features
  W <- W[with(W, order(-abs(value))), ]
  
  # Sort according to weights for each factor
  W <- as.data.frame(top_n(group_by(W, factor), n = nfeatures, wt = value))
  #
  
  # Make features names unique
  W$feature_id <- W$feature
  if ((length(unique(W$view)) > 1) && (nfeatures > 0) && (any(duplicated(W[W$factor == factors[1],]$feature_id)))) {
    message("Duplicated feature names across views, we will add the view name as a prefix")
    W$feature_id <- paste(W$view, W$feature, sep="_")
  }
  
  # In order to re-order features across multiple factors, 
  # make them unique for different factors
  W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
 

  
  W2<-merge(W,fs_met_to_merge, all.x=TRUE, by=c('feature_id' , 'view'))
  
  W2$color<-(tolower(W2$known)=='y')
  W2$color[is.na(W2$color)]=FALSE

  W2$color=factor(W2$color)
  
  p <- ggplot(W2, aes(x=.data$feature_id, y=.data$value)) +
    geom_point(size=2, aes(color=color))+
    geom_segment(aes(xend=.data$feature_id, color=color), linewidth=0.75, yend=0) +
    #scale_colour_gradient(low="grey", high="black") +
    #scale_colour_gradient(low="grey", high="black") +
    
    coord_flip() +
    labs(y="Weight") +
    
    # Theme
    theme_bw() +
    theme(
      axis.title.x = element_text(color='black'),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), hjust=1, color='black'),
      axis.text.x = element_text(color='black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      legend.position = 'none',

      legend.title = element_blank(),
      #legend.text = element_text(color="black"),
      legend.text = element_blank(),
      
      #legend.key = element_rect(fill='transparent'),
      legend.key = element_blank(),
      
      # facets
      strip.text = element_text(size=rel(1.2)),
      panel.background = element_blank(),
      panel.spacing = unit(1,"lines"),
      
      # gridlines
      panel.grid.major.y = element_blank(),
    ) +
    facet_wrap(~factor, nrow=1, scales="free")
  
  p
  if (sign=="negative") p <- p + scale_x_discrete(position = "top")
  
  # If absolute values are used, add the corresponding signs to the plot
  if (abs) {
    p <- p + 
      ylim(0,max(W$value)+0.1) + 
      geom_text(label=W$sign,y=max(W$value)+0.1, size=10)
  }
  
  return(p)
  
}

