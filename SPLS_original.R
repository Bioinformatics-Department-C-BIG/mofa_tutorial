
#BiocManager::install(c("mixOmics"))
#library(mixOmics)
#### 
# Preprocessing 


#first remember the names
library(dplyr )
library(resample)
library(mixOmics)

transpose_matrix<- function(df.aree){
  n <- df.aree$Symbol
  # transpose all but the first column (name)
  df.aree <- as.data.frame(t(df.aree[,-1]))
  colnames(df.aree) <- n
  #df.aree$myfactor <- factor(row.names(df.aree))
  return(df.aree)
  
}



preprocess_raw_data<-function(df, most_var, cut_n, ng){
  
  
  # first remove all zero entries
  if (cut_n){ df<-df[,1:cut_n] } 
  
  df<- as.data.frame(apply(df, 2, function(x) as.numeric(x)) )
  ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
  df<-df[,!ind]
  
  
  # take the most variable entries 
  if (most_var){
    n=round(dim(df)[2]/ng)
    mads<-apply(df,2,mad)
    df_selected=df[,rev(order(mads))[1:n]]
    return(df_selected)
  }else {
    return(df)
  }
  
  
  
  
  
}



####################################
#### Preliminary analysis with PCA
dir='/Users/efiathieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'

dir='E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'
output='bladder_cancer/plots/'

X1_raw<-read.csv(file = paste0(dir,'RNAseq_BladderCancer.csv' ))
X2_raw<-read.csv(file = paste0(dir,'Proteomics_BladderCancer.csv' ))
Y_raw<-read.csv(file = paste0(dir,'pheno_BladderCancer.csv' ), nrows = 16)

X1_t_raw<-transpose_matrix(X1_raw)
X2_t_raw<-transpose_matrix(X2_raw)

most_var=TRUE


ng_g=round(3,2)
ng_p=round(2,2)
X1_t<-preprocess_raw_data(X1_t_raw,most_var=most_var, cut_n=27000, ng_g)
X2_t<-preprocess_raw_data(X2_t_raw,most_var=most_var,cut_n = FALSE, ng_p)


X1_t<-log2(X1_t+1) 
X2_t<-log2(X2_t+1) 

# question: should I normalize the data? use broad institute recomendations in protigy? 


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
################
#### Preprocessing -transpose

### Preprocessing select the msot variable genes by the Median absolute deviation



ncomp_g=2
ncomp_p=5
param_str_g<-paste0( '_', most_var, '_', ncomp_g, '_ng_', round(1/ng_g,2))
param_str_g_plot<-paste( 'most var = ', most_var, ', ng =', round(1/ng_g,2),  ', ncomp = ', ncomp_g)


pca.gene <- pca(X1_t, ncomp = ncomp_g, center = TRUE, scale = TRUE)


png(paste0(output, 'pca_gene_percent', param_str_g, '.png'))
plot(pca.gene)
dev.off()



png(paste0(output,'pca_gene', param_str_g,  '.png'))
plotIndiv(pca.gene, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Bladder gene, PCA comp 1 - 2', 
          subtitle = paste0(param_str_g_plot) )

dev.off()
spca.result <- spca(X1_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 10,5))

pc_genes1<-selectVar(spca.result, comp=1)$value
pc_genes2<-selectVar(spca.result, comp=2)$value

write.csv(pc_genes1, paste0(output,'Vars_genes',param_str, '_1_X','.csv'))
write.csv(pc_genes2, paste0(output,'Vars_genes',param_str, '_2_X','.csv'))


cim(spca.result, xlab = "genes", ylab = "Samples", save='png', 
    name.save =paste0(output,'cim_genes', param_str_g), 
    title = paste0(param_str_g_plot) ) 

png(paste0(output, 'pca_genes.png'))
plotVar(spca.result, cex=c(5))
dev.off()

# Inspect row data
X1_t[,rownames(pc_genes)[1:3]]


param_str_p<-paste0( '_', most_var, '_', ncomp_p, '_ng_p', round(1/ng_p,2))
param_str_p_plot<-paste( 'most var = ', most_var, ', ng =', round(1/ng_p,2),  ', ncomp = ', ncomp_g)

ncomp_p=5
pca.proteomics <- pca(X2_t, ncomp = ncomp_p, center = TRUE, scale = TRUE)
png(paste0(output, 'pca_proteomics_percent', param_str_p, '.png'))
plot(pca.proteomics)
dev.off()
png(paste0(output, 'pca_proteomics', '_', param_str_p, '.png'))
plotIndiv(pca.proteomics, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE,
          title = paste0('Proteomics, PCA comp 1 - 2'),
          subtitle = param_str_p_plot)
dev.off()

spca.result <- spca(X2_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 10, 5))

pc_proteins1<-selectVar(spca.result, comp = 1)$value
pc_proteins2<-selectVar(spca.result, comp = 2)$value

write.csv(pc_proteins1, paste0(output,'Vars_pr',param_str_p, '_1_Y','.csv'))
write.csv(pc_proteins2, paste0(output,'Vars_pr',param_str_p, '_2_Y','.csv'))


png(paste0(output, 'pca_proteins',param_str_p, '.png'))
plotVar(spca.result, title =param_str_p_plot )
dev.off()


cim(spca.result,xlab = "proteins", ylab = "Samples",
    title = param_str_p_plot,
    save='png',
    name.save =paste0(output,'cim_proteins', param_str_p),
) 




####################################################################
##### SPLS DA







#PLS # calclate the q2 criterion used in simca-p 
### TUNING
ncomp=2
param_str<-paste0( '_', most_var, '_', ncomp_g, '_ng_g_', round(1/ng_g,2),'_ncomp_g_', ncomp_p, '_ng_p_', round(1/ng_p,2) )

param_str_plot = paste( 'most var = ', most_var, ', ng =', round(1/ng_g,2), ' ng_p = ' , round(1/ng_p,2), ', ncomp = ', ncomp_g)
result.pls <- pls(X1_t, X2_t, ncomp = ncomp, mode='canonical')  # where ncomp is the number of dimensions/components to choose
perf.pls <- perf.pls <- perf(result.pls, validation = "loo",
                             progressBar = FALSE, nrepeat = 10)

png(paste0(output,'/PLS_Q2',param_str,'.png'))
plot(perf.pls, criterion = 'Q2.total')
abline(h = 0.0975)
dev.off()


########## Tune number of features
# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 


tune.spls.bladder <- tune.spls(X1_t, X2_t, ncomp = ncomp,
                               test.keepX = list.keepX,
                               test.keepY = list.keepY,
                               nrepeat = 1, folds =4, # use 10 folds
                               mode = 'canonical', measure = 'cor')

tune.spls.bladder$choice.keepX
tune.spls.bladder$choice.keepY
# extract optimal number of variables for X dataframe
optimal$keepX <- tune.spls.bladder$choice.keepX 

#extract optimal number of variables for Y datafram
optimal$keepY <- tune.spls.bladder$choice.keepY
#optimal.keepX <- c(comp1=20, comp2=50)

#optimal.keepY <- c(comp1=10, comp2=9)

optimal$ncomp <-  length(optimal.keepX) # extract optimal number of components

fname=paste0(output,'/PLS_tune', param_str, '.png')
png(fname)
plot(tune.spls.bladder)         # use the correlation measure for tuning



fname<-paste0('bladder_cancer/settings/optimal', param_str, '.csv') 
saveRDS(optimal,fname)
readRDS(file = fname)


#result.spls <- spls(X1_t, X2_t, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
#spls.bladder<-result.spls

# use all tuned values from above
final.spls.bladder <- spls(X1_t, X2_t, ncomp = optimal$ncomp, 
                           keepX = optimal$keepX,
                           keepY = optimal$keepY,
                           mode = "canonical") # explanitory approach being used, 


plotIndiv(final.spls.bladder, ind.names = TRUE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = Y_raw$Subtype, # colour by time group
          col.per.group = color.mixo(1:2),                      # by dose group
          legend = TRUE, legend.title = 'Time', 
          legend.title.pch = 'Proteomic subtypes', 
          title = param_str_plot)



# SPLS
ncomp = 5
result.spls <- spls(X1_t, X2_t, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  criterion = 'all', progressBar = FALSE)

props<-as.character(unlist(optimal.keepY))
fname<-paste0(output,'/PLS_clustering', param_str, props, '.png')
png(fname)
plotIndiv(result.spls, ind.names = TRUE,
          rep.space = "XY-variate", # plot in Y-variate subspace,
          cex = c(7,7),
          group =  Y_raw$Subtype, # colour by time group
          legend = TRUE)
dev.off()


plotArrow(final.spls.bladder, ind.names = FALSE,
          group = Y_raw$Subtype, # colour by time group
          col.per.group = color.mixo(1:2),
          legend.title = 'Time.Group')

#########variable plots
png(paste0(output,'/PLS_correlation_circle', param_str, props, '.png'))
plotVar(final.spls.bladder , cex = c(5,5), var.names = c(TRUE, TRUE),
        title=paste0('Correlation Circle Plot',param_str), 
)
dev.off()

########network

color.edge <- color.GreenRed(50)  # set the colours of the connecting lines
common_ind<-final.spls.bladder$names$colnames$X %in%final.spls.bladder$names$colnames$Y
final.spls.bladder$names$colnames$X[common_ind]<-
  paste0(final.spls.bladder$names$colnames$X[common_ind], '__g')



#X11() # To open a new window for Rstudio
#try comp 1 or 2 
comp=1
network(final.spls.bladder, comp = comp,
        cutoff = 0.8, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = paste0(output,'sPLS Bladder Cancer Network Plot', param_str))



###### correlation plot
comp=c(1,2)
png(paste0(output,'/PLS_cim',  props, param_str, '_comp_', as.character(unlist(comp)), '.png'))
cim(final.spls.bladder, comp = comp, xlab = "proteins", ylab = "genes",
    row.cex=1.5, col.cex = 1.5, margins=c(8,8))
dev.off()


#####save important features
saveRDS(final.spls.bladder, paste0(output,'final_spls', param_str, '.RDS'))
sel_vars_1=selectVar(final.spls.bladder, comp=1)
sel_vars_2=selectVar(final.spls.bladder, comp=2)

write.csv(sel_vars_1$X, paste0(output,'Vars',param_str, '_1_X','.csv'))
write.csv(sel_vars_1$Y, paste0(output,'Vars',param_str, '_1_Y','.csv'))
write.csv(sel_vars_2$X, paste0(output,'Vars',param_str, '_2_X','.csv'))
write.csv(sel_vars_2$Y, paste0(output,'Vars',param_str, '_2_Y','.csv'))

selectVar(final.spls.bladder, comp = 2)$Y





#DIABLO
data = list(mRNA = X1_t, 
            proteomics = X2_t )


