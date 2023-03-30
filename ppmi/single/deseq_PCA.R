







# convert the columns we will be using for annotation into factors
metadata <- metadata %>%
  dplyr::mutate(
    refinebio_treatment = factor(
      refinebio_treatment,
      # specify the possible levels in the order we want them to appear
      levels = c("pre-adt", "post-adt")
    ),
    refinebio_disease = as.factor(refinebio_disease)
  )


##
all.equal(colnames(highly_variable_genes_mofa), se_filt$PATNO_EVENT_ID )

vsd$EVENT_ID
rownames(vsd)
rownames(highly_variable_genes_mofa)
vsd_filt<-vsd[rownames(vsd)%in%rownames(highly_variable_genes_mofa),]


plotPCA(
  vsd,
  intgroup = c("COHORT"), ntop=500
  
)


plotPCA(
  vsd_filt,
  intgroup = c("COHORT")
  
)

plotPCA(
  highly_variable_genes_mofa)

ggsave(paste0(outdir_s,'pca_plot.jpeg' ))


#### PCA plots

pca.data <- PCA(t(highly_variable_genes_mofa), scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))


fviz_pca_ind(pca.data)





#### 
# takes a long time, better to load datasets 
 #ddsSE <- DESeqDataSet(se_filt, 
 #                            design = ~COHORT + EVENT_ID)
 #ddsSE<-estimateSizeFactors(ddsSE)
 
 #vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
 
#rld <- rlog(ddsSE)
#plotPCA(vsd)

# also possible to perform custom transformation:
#ddsSE <- estimateSizeFactors(ddsSE)
# shifted log of normalized counts
#se <- SummarizedExperiment(vsd,
  #                         colData=colData(ddsSE))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
dim(vsd)

plotPCA( DESeqTransform( vsd ), 
         intgroup=c('COHORT', 'SEX'), ntop=20)
ggsave(paste0(outdir_s, '/pheatmap.jpeg'))






