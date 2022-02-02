install.packages('iCluster')
library("iCluster")

BiocManager::install(c("iClusterPlus"))
library(iClusterPlus)

library(GenomicRanges)
library(gplots)
library(lattice)


data
fit.single=iClusterPlus(dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics),
                        type=c( "gaussian","gaussian"),
                        lambda=c(0.04,0.61,0.90),K=1,maxiter=20)

fit.single$clusters
# Separate analysis
fit.single_genes=iClusterPlus(dt1=as.matrix(data$mRNA),
                        type=c( "gaussian"),
                        lambda=c(0.04,0.61,0.90),K=1,maxiter=20)

fit.single_genes$clusters

fit.single_proteomics=iClusterPlus(dt1=as.matrix(data$proteomics),
                        type=c( "gaussian"),
                        lambda=c(0.04,0.61,0.90),K=1,maxiter=20)

fit.single_proteomics$clusters

#### MODEL TUNING - k - number of clusters
for(k in 1:2){
  cv.fit = tune.iClusterPlus(cpus=1,dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics), 
                    type=c("gaussian","gaussian"),K=k, n.lambda=8,
      scale.lambda=c(1,1,1),maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
  }



######## model selection 

output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

