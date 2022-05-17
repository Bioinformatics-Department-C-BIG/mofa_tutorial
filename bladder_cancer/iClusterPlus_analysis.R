
source('bladder_cancer/iClusterPlus.R')




cv.fit$lambda

# modes: prot, gene, multi 
mode='prot'

######## model selection 
n.lambda = 217
output=alist()
if (mode=='prot'){
  files=grep(paste0("cv.fit.prot", n.lambda),dir(settings))
  
}
if (mode=='multi'){
  files=grep(paste0("cv.fit.k_new_de_data_cluster", n.lambda),dir(settings))
  
  }



for(i in 1:length(files)){
  load(paste0(settings, dir(settings)[files[i]]))
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output) ; nK
BIC = getBIC(output); BIC
devR = getDevR(output) ; devR

minBICid = apply(BIC,2,which.min); minBICid
devRatMinBIC = rep(NA,nK)

for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}



# The optimal k (number of latent variables) is where the curve of %Explained variation levels off. 
# Examine this to select the optimal k 
png(paste0(icluster_out, paste0('iCluster_explained_var', n.lambda, mode, '.png')))
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
dev.off()
# Proteomics k=2
k=2


clusters=getClusters(output)
rownames(clusters)=rownames(data$mRNA)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]] # best cluster is the one with min BIC
getwd()
icluster_out<-'bladder_cancer/iCluster/'





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
for(i in 1:length(best.fit$beta)){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("expression", "proteomics")


#print a few examples of selected features
write.csv(sigfeatures[[1]], paste0(icluster_out,'Vars_genes',param_str_icl, '_X','.csv'))
write.csv(sigfeatures[[2]], paste0(icluster_out,'Vars_proteins',param_str_icl, '_X','.csv'))


head(sigfeatures[[1]], 15)
sigfeatures[[1]]
sigfeatures[[1]] %in% colnames(data$mRNA)
  head(sigfeatures[[2]])


## PROBLEM HERE

bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bluered(256)
col.scheme[[2]] = bluered(256)


if (mode == 'gene'){
  png(paste0(icluster_out,param_str_icl,'_k_', K_genes, 'heatmap_genes.png'))
  plotHeatmap(fit=fit.single_genes,datasets=list(mRNA.image),
              type=c("gaussian","gaussian"), col.scheme = col.scheme,
              row.order=c(T,T),sparse=c(T,T),
              cap=c(F,F))
  dev.off()
}

if (mode=='prot'){
  fname = paste0(icluster_out,'Vars_genes',param_str_icl, '_X', n.lambda,'.csv')
  write.csv(sigfeatures[[1]], fname)
  
  png(paste0(icluster_out,param_str_icl,'_k_', K_proteomics, 'heatmap_proteins.png'))
  plotHeatmap(fit=best.fit,datasets=list(prot.image),
              type=c("gaussian","gaussian"), col.scheme = col.scheme,
              row.order=c(T,T),sparse=c(T,T),
              cap=c(F,F))
  dev.off()
  
}



if (mode=='multi'){
  
png(paste0(icluster_out,param_str_icl, 'heatmap_multi_trained.png'))
plotHeatmap(fit=best.fit,datasets=list(mRNA.image,prot.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()


} 


png(paste0(icluster_out,param_str_icl, 'heatmap_multi.png'))

plotHeatmap(fit=fit.single,datasets=list(mRNA.image,prot.image),
            type=c("gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(T,T,T),sparse=c(T,T),
            cap=c(F,F))
dev.off()



### TODO Need to reanaluze with new fitting parameters 


