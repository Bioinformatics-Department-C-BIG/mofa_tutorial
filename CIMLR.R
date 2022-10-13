# run this commands only if the following R packages are not already installed
install.packages("devtools", dependencies = TRUE)
install.packages("Matrix", dependencies = TRUE)


# install CIMLR from Github
library("devtools")
install_github("danro9685/CIMLR", ref = 'R')

# load CIMLR library
library("CIMLR")



library(CIMLR)
data(GliomasReduced)


set.seed(11111)
NUMC = 2:10
res_example = CIMLR_Estimate_Number_of_Clusters(GliomasReduced$in_X,
                                                NUMC = NUMC,
                                                cores.ratio = 0)


set.seed(11111)
example = CIMLR(X = GliomasReduced$in_X, c = 5, cores.ratio = 0)


plot(example$ydata,
     col = c(topo.colors(5))[example$y[["cluster"]]],
     xlab = "CIMLR component 1",
     ylab = "CIMLR component 2",
     pch = 20,
     main="CIMILR 2D visualization for GliomasReduced")


set.seed(11111)
input_data = rbind(GliomasReduced$in_X$point_mutations,GliomasReduced$in_X$copy_numbers,
                   GliomasReduced$in_X$methylations,GliomasReduced$in_X$expression_values)
ranks = CIMLR_Feature_Ranking(A=example$S,X=input_data)


head(ranks$pval)

head(ranks$aggR)



#### bladder cancer 


data_in_cimlr = list(mRNA = highly_variable_genes_voom, 
            proteomics = highly_variable_proteins_voom )

example = CIMLR(X =data_in_cimlr , c = 3, cores.ratio = 0)

input_data = rbind(highly_variable_genes_voom,highly_variable_proteins_voom)
ranks = CIMLR_Feature_Ranking(A=example$S,X=input_data)


head(ranks$pval)

head(ranks$aggR,80)


example$y[["cluster"]]


plot(example$ydata,
     col = c(topo.colors(5))[example$y[["cluster"]]],
     xlab = "CIMLR component 1",
     ylab = "CIMLR component 2",
     pch = 20,
     main="CIMILR 2D visualization for GliomasReduced")



### Check which features are most important ? 
rownames(input_data)[head(ranks$aggR,80)]

