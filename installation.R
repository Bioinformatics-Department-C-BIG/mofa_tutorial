#install.packages("BiocManager", repos = "http://cran.us.r-project.org")
install.packages('devtools')


BiocManager::install("mixOmics")
BiocManager::install(version = "3.16")
#remove.packages('BiocManager')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version(16))
remove.packages('BiocManager')
BiocManager::install("MOFA2")

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

#BiocManager::install('MultiAssayExperiment')
install.packages("remotes")
remotes::install_github("bioFAM/MOFA2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


#BiocManager::install("MOFAdata")'

#remove.packages('MOFA2')

#To use the latest features of MOFA you can install the software from GitHub:
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
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
#install.packages('ggnewscale')
BiocManager::install('GOfuncR')
BiocManager::install('ensembldb')
BiocManager::install('org.Hs.eg.db')




install.packages("remotes")
remotes::install_github("jmw86069/jamenrich")

