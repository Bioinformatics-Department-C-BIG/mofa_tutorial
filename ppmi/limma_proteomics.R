cran_packages <- c("remotes", "dplyr", "ggplot2")
for (pkg_i in cran_packages) {
  if (!require(pkg_i, quietly = T, character.only = T))
    install.packages(pkg_i)
}


remotes::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
## ------------------------
library(MSnSet.utils)
library(dplyr)
library(ggplot2)

# MSnSet for testing
data("cptac_oca")
m <- oca.set
mnew<-m@assayData$exprs
