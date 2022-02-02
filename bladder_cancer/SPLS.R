
BiocManager::install(c("mixOmics"))
library(mixOmics)
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



preprocess_raw_data<-function(df, most_var, cut_n){
      
      
      # first remove all zero entries
      if (cut_n){ df<-df[,1:cut_n] } 
  
      df<- as.data.frame(apply(df, 2, function(x) as.numeric(x)) )
      ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
      df<-df[,!ind]
      
      
      # take the most variable entries 
      if (most_var){
        n=round(dim(df)[2]/4)
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
X1_t<-preprocess_raw_data(X1_t_raw,most_var=most_var, cut_n=27000)
X2_t<-preprocess_raw_data(X2_t_raw,most_var=most_var,cut_n = FALSE)


X1_t<-log2(X1_t+1) 
X2_t<-log2(X2_t+1) 

# question: should I normalize the data? use broad institute recomendations in protigy? 


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
################
#### Preprocessing -transpose

### Preprocessing select the msot variable genes by the Median absolute deviation




ncomp_g=2
pca.gene <- pca(X1_t, ncomp = ncomp_g, center = TRUE, scale = TRUE)
png(paste0(output, 'pca_gene_percent', '_', most_var, '_', ncomp_g, '.png'))
plot(pca.gene)
dev.off()


png(paste0(output,'pca_gene', '_', most_var, '.png'))
plotIndiv(pca.gene, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Bladder gene, PCA comp 1 - 2')

dev.off()
spca.result <- spca(X1_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 10,5))

pc_genes<-selectVar(spca.result, comp = 1)$value
cim(spca.result, xlab = "genes", ylab = "Samples", save='png', 
    name.save =paste0(output,'cim_genes') ) 

png(paste0(output, 'pca_genes.png'))
plotVar(spca.result, cex=c(5))
dev.off()

# Inspect row data
X1_t[,rownames(pc_genes)[1:3]]


ncomp_p=5
pca.proteomics <- pca(X2_t, ncomp = ncomp_p, center = TRUE, scale = TRUE)
png(paste0(output, 'pca_proteomics_percent', '_', most_var,'_', ncomp_p, '.png'))
plot(pca.proteomics)
dev.off()
png(paste0(output, 'pca_proteomics', '_', most_var, '.png'))
plotIndiv(pca.proteomics, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE,
          title = paste0('Proteomics, PCA comp 1 - 2'),
          subtitle = paste0( 'most_var = ', most_var))
dev.off()

spca.result <- spca(X2_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 5, 5))

selectVar(spca.result, comp = 1)$value

png(paste0(output, 'pca_proteins.png'))
plotVar(spca.result)
dev.off()

cim(spca.result, xlab = "proteins", ylab = "Samples",
    save='png',name.save =paste0(output,'cim_proteins') ) 




####################################################################
##### SPLS DA







#PLS # calclate the q2 criterion used in simca-p 
### TUNING

ncomp=2
result.pls <- pls(X1_t, X2_t, ncomp = ncomp, mode='canonical')  # where ncomp is the number of dimensions/components to choose
perf.pls <- perf.pls <- perf(result.pls, validation = "loo",
                             progressBar = FALSE, nrepeat = 10)

png(paste0(output,'/PLS_Q2.png'))
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

png(paste0(output,'/PLS_tune.png'))
tune.spls.bladder$choice.keepX
tune.spls.bladder$choice.keepY
# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.bladder$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.bladder$choice.keepY
optimal.keepY <- c(comp1=10, comp2=9)

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components


plot(tune.spls.bladder)         # use the correlation measure for tuning
dev.off()



#result.spls <- spls(X1_t, X2_t, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
#spls.bladder<-result.spls

# use all tuned values from above
final.spls.bladder <- spls(X1_t, X2_t, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "canonical") # explanitory approach being used, 


plotIndiv(final.spls.bladder, ind.names = TRUE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = Y_raw$Subtype, # colour by time group
          col.per.group = color.mixo(1:2),                      # by dose group
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Proteomic subtypes')



# SPLS
ncomp = 5
result.spls <- spls(X1_t, X2_t, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  criterion = 'all', progressBar = FALSE)

props<-as.character(unlist(optimal.keepY))
png(paste0(output,'/PLS_clustering', props, '.png'))
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
png(paste0(output,'/PLS_correlation_circle',props, '.png'))
plotVar(final.spls.bladder , cex = c(4,4), var.names = c(TRUE, TRUE))
dev.off()

########network

color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

common_ind<-final.spls.bladder$names$colnames$X %in%final.spls.bladder$names$colnames$Y
final.spls.bladder$names$colnames$X[common_ind]<-
  paste0(final.spls.bladder$names$colnames$X[common_ind], '__g')



#X11() # To open a new window for Rstudio
network(final.spls.bladder, comp = 1:2,
        cutoff = 0.8, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = paste0(output,'sPLS Bladder Cancer Network Plot'))



###### corelation plot
png(paste0(output,'/PLS_cim',  'props', '.png'))
cim(final.spls.bladder, comp = 1:2, xlab = "proteins", ylab = "genes",
    row.cex=1.5, col.cex = 1.5, margins=c(8,8))
dev.off()


selectVar(final.spls.bladder, comp = 1)$Y
selectVar(final.spls.bladder, comp = 2)$Y


#### TUNING 
# Enrichment analysis

reticulate::use_python(python = "C:/Users/athienitie/Anaconda3/python.exe")


library(mixOmics)
data(breast.TCGA)
#BiocManager::install("biomaRt")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

#### Enrichment analysis 
utils::data()



# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)



#DIABLO
data = list(mRNA = X1_t, 
            proteomics = X2_t )

Y<-Y_raw$Subtype


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

design 

# tune components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                         design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10)

#perf.diablo  # lists the different outputs
plot(perf.diablo) 


sgccda.res$design
selectVar(sgccda.res, block = 'mRNA', comp = 1)$mRNA$name 


plotDiablo(sgccda.res, ncomp = 1)
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')



plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

dev.off()  
plotVar(sgccda.res, cutoff = 0.7, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'brown1' ))


circosPlot(sgccda.res, cutoff = 0.8,line = TRUE, 
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3", ), size.labels = 1.5)

