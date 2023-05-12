# install.packages('iCluster')
library("iCluster")
library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
# BiocManager::install('iClusterPlus')
tune.iClusterBayes(cpus=6,dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics),
type=c("gaussian","gaussian"),
K=1:6)
library('iClusterPlus')
tune.iClusterBayes(cpus=6,dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics),
type=c("gaussian","gaussian"),
K=1:6)

icluster_results<-tune.iClusterBayes(cpus=6,dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics),
K=1:6)
icluster_results
icluster_results
icluster_results$fit[[1]]


# todo: add the clusters onto the analysis
