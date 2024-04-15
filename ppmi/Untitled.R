
install.packages('minfi')
BiocManager::install('minfiData')

library(minfi)
library(minfiData)
## RGsetEx: RGChannelSet, 622,399 features
MsetEx <- preprocessRaw(RGsetEx)
## MsetEx: MethylSet, 485,512 features
GMsetEx <- mapToGenome(MsetEx)
## GMsetEx: GenomicMethylSet, 485,512 features

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)


methylation_files<-read.csv('ppmi/data/ppmi_140_link_list_20210607.csv')
unique_patients<-unique(methylation_files$PATNO)
## how many instances for each aptients??

