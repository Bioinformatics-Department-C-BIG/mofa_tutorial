
BiocManager::install('compGenomRData')
install.packages('devtools')
devtools::install_github("compgenomr/compGenomRData")
devtools::install_github("compgenomr/compGenomRData")

### tutorial https://compgenomr.github.io/book/use-case-multi-omics-data-from-colorectal-cancer.html

library(compGenomRData)
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")



library(FactoMineR)


r.mfa <- FactoMineR::MFA(
  t(rbind(x1,x2,x3)), # binding the omics types together
  c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
  graph=FALSE)
Since this generates a two-dimensional factorization of the multi-omics data, we can now plot each tumor as a dot in a 2D scatter plot to see how well the MFA factors separate the cancer subtypes. The following code snippet generates Figure 11.6:
  
  
  # first, extract the H and W matrices from the MFA run result
  mfa.h <- r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord

# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- factor(covariates[rownames(mfa_df),]$cms_label)

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
  ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")

