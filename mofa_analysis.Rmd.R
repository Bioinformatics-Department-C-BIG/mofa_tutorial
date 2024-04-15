

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained','.png'), width = 4, height=4, dpi=100)


MOFAobject@samples_metadata$Grade<-as.factor(MOFAobject@samples_metadata$Grade)
MOFAobject@samples_metadata$Grade
MOFAobject@samples_metadata
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("Grade", "EORTC.risk", 'Subtype'), 
                                  plot = "log_pval",
                                  
)
ggsave(paste0(outdir, 'factors_covariates','.png'), width = 4, height=4, dpi=100)


TOP<-MOFA2::get_weights(MOFAobject,views = 'mRNA', factors=1, 
                        as.data.frame =TRUE)


TOP[1:100,]
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


#install.packages('psych')

plot_factors(MOFAobject, 
             factors = c(1,2), 
             dot_size = 2.5
)


plot_weights(MOFAobject,
             view = "proteomics",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)




### Age, gender , stage does not discriminate factors
#Conclusion here:# factor 1 correlates with grade

p<-plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "EORTC.risk",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)
p
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Subtype",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Subtype",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)


plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Age",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

##### plot weights 
fps<-seq(1:5)

vps<-c("proteomics" ,"mRNA")

v_set=c()
v_set=c()
fps<-seq(1:3)

plot_top_weights(MOFAobject,
                 view = 'mRNA',
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave(paste0(outdir, 'top_weights_','2', 'mRNA','_','.png'), width = , height=4, dpi=100)



for (i in 1:length(vps)){
  for (ii in 1:length(fps)){
    
    plot_top_weights(MOFAobject,
                     view = vps[i],
                     factor = fps[ii],
                     nfeatures = 10,     # Top number of features to highlight
                     scale = T           # Scale weights from -1 to 1
    )
    ggsave(paste0(outdir, 'top_weights_', fps[ii],'_',vps[i],'.png'), width = , height=4, dpi=100)
    
    plot_weights(MOFAobject, 
                 view = vps[i], 
                 factor = fps[ii], 
                 nfeatures = 10
    )
    ggsave(paste0(outdir, 'all_weights_', fps[ii],'_',vps[i],'.png'), width = 4, height=4, dpi=100)
    
    
    
    ###### Heatmaps 
    nfs=20
    
    jpeg(paste0(outdir, 'heatmap_',fps[ii],'_','mRNA_', 'nfs_', nfs, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
    fps[ii]=1
    # Plot heatmaps for each factor only for mRNA 
    p<-plot_data_heatmap(MOFAobject, 
                         view = vps[i], 
                         factor =  fps[ii],  
                         features = nfs,
                         denoise = TRUE,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         scale = "row",
                         fontsize_number = 5
                         
                         
    )
    dev.off()
    
    jpeg(paste0(outdir, 'heatmap_',fps[ii],'_','Drugs_', 'nfs_', nfs, '.jpeg'), height=20*nfs, width=20*nfs)
    
    p<-plot_data_heatmap(MOFAobject, 
                         view = vps[i], 
                         factor =  fps[ii],  
                         features = nfs,
                         denoise = TRUE,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         scale = "row"
    )
    dev.off()
    
    # top weights
    # concat all 
    
    
    
  }
  
  
}








plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


### So multi omics factors are more related to Stage than to subtype!!
##How?? plot on the factors and color by stage

p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "EORTC.risk",
                  shape_by = "Grade",
                  dot_size = 2.5,
                  show_missing = T
)

#p <- p + 
#  geom_hline(yintercept=-1, linetype="dashed") +
#  geom_vline(xintercept=(-0.5), linetype="dashed")
print(p)
ggsave(p)
ggsave(paste0(outdir,'factor_plot','.png'), width = 4, height=4, dpi=120)



### 
## GENE SET ENRICHMENT! 
## AND reactome gsa enrichment!!


###### GSEA 
library(MOFAdata)

utils::data(reactomeGS)

head((reactomeGS))

features_names(MOFAobject)$mRNA
# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)



# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "negative"
)

for (ii in 1:4){
  sign<-res.positive$pval.adj[,ii][which(res.positive$pval.adj[,ii]<0.05)]
  #sign<-res.negative$pval.adj[,ii][which(res.negative$pval.adj[,ii]<0.005)]
  
  print(length(sign))
}
names(res.positive)[which(res.positive$pval.adj<0.05)]

View(res.positive$pval.adj)


theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

names(res.positive)
dev.off()
# ENS ids seems to be required for the enrichment 
plot_enrichment_heatmap(res.positive)
ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)
dev.off()

plot_enrichment_heatmap(res.negative)
ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


# Positive Factor 1: Hypoxia, oxygen dependent etc.
names(which(res.positive$pval.adj[,'Factor1']<0.0005))[1:10]
names(which(res.negative$pval.adj[,'Factor1']<0.05))

# Positive factor 2: citric acid cycle, meabolism, L13a-mediated, respratory
names(which(res.positive$pval.adj[,'Factor2']<0.00000005))[1:10]
names(which(res.negative$pval.adj[,'Factor2']<0.05))


# Positive factor 3: chromatin modifying enzymes, regulation of tp53
names(which(res.positive$pval.adj[,'Factor3']<0.05))[1:10]
names(which(res.negative$pval.adj[,'Factor3']<0.05))

#  
names(which(res.positive$pval.adj[,'Factor4']<0.05))[1:5]
names(which(res.negative$pval.adj[,'Factor4']<0.05))[1:5]

