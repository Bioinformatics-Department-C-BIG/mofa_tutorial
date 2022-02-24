#install.packages('iCluster')
library("iCluster")

#BiocManager::install(c("iClusterPlus"))
library(iClusterPlus)

library(GenomicRanges)
library(gplots)
library(lattice)


#DIABLO
data = list(mRNA = X1_t , 
            proteomics = X2_t )


##NOTE ON pre processing: LOG first then standardize!!                                                  
data = list(mRNA = t(highly_variable_genes_voom), 
            proteomics = t(highly_variable_proteins_voom) )


param_str_icl<-paste0( 'DE_data','_', most_var, '_ng_g_', round(1/ng_g,2),'_ncomp_g_','_ng_p_', round(1/ng_p,2)  )



# Separate analysis
K_genes=2
fit.single_genes=iClusterPlus(dt1=as.matrix(data$mRNA),
                              type=c( "gaussian"),
                              lambda=c(0.04,0.61,0.90),K=K_genes,maxiter=20)

fit.single_genes$clusters


K_proteomics=2
fit.single_proteomics=iClusterPlus(dt1=as.matrix(data$proteomics),
                                   type=c( "gaussian"),
                                   lambda=c(0.04,0.61,0.90),K=K_proteomics,maxiter=20)

fit.single_proteomics$clusters



fit.single=iClusterPlus(dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics),
                          type=c("gaussian","gaussian"),
                          lambda=c(0.04,0.61,0.90),K=2,maxiter=10)

settings<-'bladder_cancer/settings/'
#### MODEL TUNING - k - number of clusters
## And save these data
n.lambda<-89
for(k in 1:3){
  cv.fit = tune.iClusterPlus(cpus=1,dt1=as.matrix(data$mRNA),dt2=as.matrix(data$proteomics), 
                             type=c("gaussian","gaussian"),K=k, n.lambda=n.lambda,
                             scale.lambda=c(1,1,1),maxiter=20)
  save(cv.fit, file=paste(settings,"cv.fit.k_new_de_data",n.lambda,k,".Rdata",sep=""))
}

cv.fit$lambda

######## model selection 

output=alist()
files=grep("cv.fit.k_new_de",dir(settings))
for(i in 1:length(files)){
  load(paste0(settings, dir(settings)[files[i]]))
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

clusters=getClusters(output)
rownames(clusters)=rownames(data$mRNA)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]
getwd()
icluster_out<-'bladder_cancer/iCluster/'
png(paste0(icluster_out, 'iCluster_explained_var.png'))
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
dev.off()


#### Generate heatmap

#truncate the values for a better image plot
mRNA.image=data$mRNA
mRNA.image[mRNA.image>10]=10
mRNA.image[mRNA.image< -1]= -1
prot.image=data$proteomics
prot.image[prot.image>10]=10
prot.image[prot.image< -1]= -1

mRNA.image=data$mRNA
head(summary(head(mRNA.image)))
#mRNA.image[mRNA.image>10]=10
#mRNA.image[mRNA.image< -1]= -1
prot.image=data$proteomics
#prot.image[prot.image>10]=10
#prot.image[prot.image< -1]= -1


###### Feature Selection


features = alist()
features[[1]] = colnames(data$mRNA)
features[[2]] = colnames(data$proteomics)
sigfeatures=alist()
for(i in 1:2){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("expression", "proteomics")


#print a few examples of selected features
write.csv(namessigfeatures[[1]], paste0(icluster_out,'Vars_genes',param_str_icl, '_X','.csv'))
write.csv(sigfeatures[[2]], paste0(icluster_out,'Vars_proteins',param_str_icl, '_X','.csv'))


head(sigfeatures[[1]])
sigfeatures[[1]] %in% colnames(data$mRNA)
head(sigfeatures[[2]])


## PROBLEM HERE

bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bluered(256)
col.scheme[[2]] = bluered(256)



png(paste0(icluster_out,param_str_icl,'_k_', K_genes, 'heatmap_genes.png'))
plotHeatmap(fit=fit.single_genes,datasets=list(mRNA.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()

png(paste0(icluster_out,param_str_icl,'_k_', K_proteomics, 'heatmap_proteins.png'))
plotHeatmap(fit=fit.single_proteomics,datasets=list(prot.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()


png(paste0(icluster_out,param_str_icl, 'heatmap_multi_trained.png'))
plotHeatmap(fit=best.fit,datasets=list(mRNA.image,prot.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()
png(paste0(icluster_out,param_str_icl, 'heatmap_multi.png'))

plotHeatmap(fit=fit.single,datasets=list(mRNA.image,prot.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()



### TODO Need to reanaluze with new fitting parameters 


