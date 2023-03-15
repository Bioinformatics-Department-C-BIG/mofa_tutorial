
PVS<-read.csv('/Users/efiathieniti/Downloads/covariate_corelations.csv', row.names = 1)


PVS
sig<-apply(PVS,2, function(x) {x>abs(log10(0.05/68))})

round(PVS[1,][which(sig[1,])], digits = 2)

ll<-  apply(sig,1,function(x) { x[which(x)]})
names(unlist(ll))
