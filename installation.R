install.packages("BiocManager", repos = "http://cran.us.r-project.org")
install.packages('devtools')


BiocManager::install("mixOmics")
BiocManager::install(version = "3.16")
#remove.packages('BiocManager')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version(16))
remove.packages('BiocManager')
BiocManager::install("MOFA2")
BiocManager::install("MOFAdata")

BiocManager::install('Org.Mm.eg.db')
#remove.packages('MOFA2')


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


##setup
install.packages('rstudioapi')


###### Single layer



install.packages('R.filesets')
BiocManager::install('DESeq2')
BiocManager::install("SummarizedExperiment")
install.packages('data.table')
install.packages('dplyr')
BiocManager::install('EnhancedVolcano')

install.packages("factoextra")





install.packages('tidyverse')
install.packages('UpSetR')

#### MOFA
BiocManager::install('MOFAdata')
BiocManager::install('MOFA2')



#BiocManager::install('MultiAssayExperiment')
install.packages("remotes")
remotes::install_github("bioFAM/MOFA2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("MOFAdata")
BiocManager::install("MOFA")

#remove.packages('MOFA2')

#To use the latest features of MOFA you can install the software from GitHub:
#devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
#If you do so, you have to manually install the Python dependencies using pip (from the Unix terminal). Importantly, this has to be done before the R installation.

#pip install mofapy2
### proteins
BiocManager::install('DEP')
BiocManager::install('vsn')
BiocManager::install('MultiAssayExperiment')
BiocManager::install('rbioapi')



### Enrichment 

BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
install.packages('europepmc')
BiocManager::install('tidyverse')


BiocManager::install('GOfuncR')
BiocManager::install('ensembldb')
BiocManager::install('org.Hs.eg.db')



### MIXOMICS

BiocManager::install(c("mixOmics"))
install.packages('visNetwork')
install.packages('igraph')




install.packages("remotes")
remotes::install_github("jmw86069/jamenrich")
remotes::install_github("jmw86069/jamenrich/multienrichjam")

#remotes::install_github("jmw86069/jamenrich")






#### 


install.packages('NMF')



BiocManager::install('WGCNA')
BiocManager::install('OmnipathR', force=TRUE)
BiocManager::install('dnet')
BiocManager::install('sgof')
BiocManager::install('factoextra')
BiocManager::install('GOfuncR')

BiocManager::install('VennDiagram')

BiocManager::install('reshape')

BiocManager::install('randomForest')
BiocManager::install('DescTools')

install.packages('caret')
install.packages('tibble')
suppressPackageStartupMessages(library(randomForest))



BiocManager::install("vsn")
BiocManager::install('GGally')



## install timeOmics
BiocManager::install('timeOmics')
BiocManager::install('propr')



