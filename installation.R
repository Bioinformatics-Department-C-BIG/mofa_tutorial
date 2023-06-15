install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install("mixOmics")
install.packages("BiocManager")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")

BiocManager::install('Org.Mm.eg.db')



install.packages(
  "ggplot2",
  repos = c("http://rstudio.org/_packages",
            "http://cran.rstudio.com")
)
install.packages(
  "dplyr",
  repos = c("http://rstudio.org/_packages",
            "http://cran.rstudio.com")
)

install.packages('iCluster')
library("iCluster")

BiocManager::install(c("iClusterPlus"))
BiocManager::install(version = '3.15')
BiocManager::install("iClusterPlus")


BiocManager::install("GenomicRanges")
BiocManager::install("Glimma")


install.packages('gplots')
install.packages('lattice')


install.packages('iCluster')


BiocManager::install('edgeR')
BiocManager::install('cluster')




M###### Single layer



install.packages('R.filesets')
BiocManager::install('DESeq2')
BiocManager::install("SummarizedExperiment")
install.packages('data.table')
install.packages('dplyr')
BiocManager::install('EnhancedVolcano')



install.packages('UpSetR')

#### MOFA
BiocManager::install('MOFAdata')
<<<<<<< HEAD
BiocManager::install('MOFA2')

=======

#BiocManager::install('MultiAssayExperiment')
#BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")'

### proteins
BiocManager::install('DEP')
BiocManager::install('vsn')
>>>>>>> 6e75bc7543129daa3d6c96166d5b27e2d9ac433b
BiocManager::install('MultiAssayExperiment')



### Enrichment 

BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
#install.packages('ggnewscale')
BiocManager::install('GOfuncR')
BiocManager::install('ensembldb')
BiocManager::install('org.Hs.eg.db')



### MIXOMICS

BiocManager::install(c("mixOmics"))

