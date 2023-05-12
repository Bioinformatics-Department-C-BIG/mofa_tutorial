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




###### Single layer



install.packages('R.filesets')
BiocManager::install('DESeq2')
BiocManager::install("SummarizedExperiment")
install.packages('data.table')
install.packages('dplyr')


install.packages('UpSetR')

#### MOFA
BiocManager::install('MOFAdata')




### Enrichment 

BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
#install.packages('ggnewscale')
BiocManager::install('GOfuncR')
BiocManager::install('ensembldb')
BiocManager::install('org.Hs.eg.db')

