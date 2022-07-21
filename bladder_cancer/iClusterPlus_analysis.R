
source('bladder_cancer/iClusterPlus.R')
icluster_out='bladder_cancer/iCluster/' 



cv.fit$lambda

# modes: prot, gene, multi 
mode='gene'

######## model selection 
files=list()

mode='protein'
n.lambda=13
ng_p=round(5,2)
ng_g=round(35,2)


mode='gene'
n.lambda=89
ng_p=round(5,2)
ng_g=round(35,2)


n.lambda_list=c(144,34,89,233,377,89, 217,13)
#n.lambda_list=c(13,34,89,217,89)

mode='multi'
n.lambda=233
ng_p=round(4,2)
ng_g=round(25,2)
#param_str_icl<-paste0( 'DE_data','_', most_var, '_ng_g_', round(1/ng_g,2),'_ncomp_g_','_ng_p_', round(1/ng_p,2)  )
param_str_icl<-paste0(  most_var, '_ng_g_', ng_g,'_ng_p_', ng_p  )



for (ii in 1:4){
  print(n.lambda_list[ii])
  get_results_for_run(n.lambda_list[ii])
  dev.off()
}



get_results_for_run<-function(n.lambda){
  files=list()
  output=alist()
  if (mode=='prot'){
    files=grep(paste0("cv.fit.prot", n.lambda),dir(settings))
    
  }
  if (mode=='gene'){
    files=grep(paste0("cv.fit.gene", n.lambda),dir(settings))
    files=grep(paste0("cv.fit.gene", n.lambda,'ng_', ng_g ),dir(settings))
    
  }
  if (mode=='multi'){
    files=grep(paste0("cv.fit.multi", n.lambda, 'ng_', ng_g, '_', ng_p),dir(settings))
    
    #files=grep(paste0("cv.fit.multi", n.lambda),dir(settings))
    
    }
  
  
  
  for(i in 1:length(files)){
    print(paste0(settings, dir(settings)[files[i]]))
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
  print(devRatMinBIC, digits = 4)
  
  
  # The optimal k (number of latent variables) is where the curve of %Explained variation levels off. 
  # Examine this to select the optimal k 
  png(paste0(icluster_out, paste0('iCluster_explained_var', param_str_icl, n.lambda, mode, '.png')))
  plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
       ylab="%Explained Variation")
  dev.off()
  # Proteomics k=2
  
  
  
  
  clusters=getClusters(output)
  rownames(clusters)=rownames(data$mRNA)
  colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
  #write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
  k=which.max(devRatMinBIC)
  k=2
  for (k in 2:nK){
    best.cluster=clusters[,k]
    best.fit=output[[k]]$fit[[which.min(BIC[,k])]] 
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
    #head(summary(head(mRNA.image)))
    #mRNA.image[mRNA.image>10]=10
    #mRNA.image[mRNA.image< -1]= -1
    prot.image=data$proteomics
    #prot.image[prot.image>10]=10
    #prot.image[prot.image< -1]= -1
    
    
    ###### Feature Selection
    sigfeatures=alist()
    features=alist()
    features[[1]] = colnames(data$mRNA)
    features[[2]] = colnames(data$proteomics)
    
    if (mode=='gene'){
      features[[1]] = colnames(data$mRNA)
      
    }else if (mode=='prot'){
      features[[1]] = colnames(data$proteomics)
      
    }
    
    upper_quantile=0.9
    for(i in 1:length(best.fit$beta)){
      rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
      upper=quantile(rowsum,prob=upper_quantile)
      sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
      
    }
    
    
    #print a few examples of selected features
    
  
    if (length(sigfeatures)>1){
      names(sigfeatures)=c("expression", "proteomics")
      
      write.csv(sigfeatures[[2]], paste0(icluster_out,'Vars_proteins',param_str_icl, '_X','.csv'))
      
    }
  
    
    head(sigfeatures[[1]], 15)
    sigfeatures[[1]]
    sigfeatures[[1]] %in% colnames(data$mRNA)
  
    
    ## PROBLEM HERE
    
    bw.col = colorpanel(2,low="white",high="black")
    col.scheme = alist()
    col.scheme[[1]] = bluered(256)
    col.scheme[[2]] = bluered(256)
    
    
    if (mode == 'gene'){
      fname = paste0(icluster_out,'Vars_genes',param_str_icl, '_X_', k ,n.lambda)
      write.csv(sigfeatures[[1]], paste0(fname, '_', upper_quantile,'.csv' ))
      top_genes=sigfeatures[[1]]
      top_genes
      gene_clusters=best.cluster
      write.csv(best.cluster, paste0(fname,  'best_clusters.csv'))
      
      png(paste0(icluster_out,param_str_icl,'_k_', k, '_' ,n.lambda, '_heatmap_genes.png'))
      plotHeatmap(fit=best.fit,datasets=list(mRNA.image),
                  type=c("gaussian","gaussian"), col.scheme = col.scheme,
                  row.order=c(T,T),sparse=c(T,T),
                  cap=c(F,F))
      dev.off()
    }
    
    if (mode=='prot'){
      fname = paste0(icluster_out,'Vars_prot',param_str_icl, '_X_', k, n.lambda)
      write.csv(sigfeatures[[1]], paste0(fname, '_', upper_quantile, '.csv' ))
      write.csv(best.cluster, paste0(fname,  'best_clusters.csv'))
      top_proteins=sigfeatures[[1]]
      protein_clusters=best.cluster
      
      
      png(paste0(icluster_out,param_str_icl,'_k_', k, '_' , n.lambda, '_heatmap_proteins.png'))
      plotHeatmap(fit=best.fit,datasets=list(prot.image),
                  type=c("gaussian","gaussian"), col.scheme = col.scheme,
                  row.order=c(T,T),sparse=c(T,T),
                  cap=c(F,F))
      dev.off()
      
    }
    
    
    
    if (mode=='multi'){
    fname = paste0(icluster_out,'Vars_multi',param_str_icl, '_X_', k,  n.lambda)
    write.csv(sigfeatures[[1]], paste0(fname,'_', upper_quantile, '.csv' ))
    write.csv(best.cluster, paste0(fname,  'best_clusters.csv'))
    multi_clusters=best.cluster
    
    
   # png(paste0(icluster_out,param_str_icl, 'heatmap_multi_trained.png'))
   # plotHeatmap(fit=best.fit,datasets=list(mRNA.image,prot.image),
  #              type=c("gaussian","gaussian"), col.scheme = col.scheme,
  #              row.order=c(T,T),sparse=c(T,T),
  #              cap=c(F,F))
  #  dev.off()
    
    
    } 
    
  }
  
  
  


}


ha<- heatmap(mRNA.image)
  
  

### TODO Need to reanalyze with new fitting parameters 


