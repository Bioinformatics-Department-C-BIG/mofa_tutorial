head(normtissue <- hpaNormalTissue())
colnames(normtissue)
ens_orig1
normtissue$Gene.name


normtissue_matched<-normtissue[match(ens_orig1, normtissue$Gene ), ]

normtissue_matched
ens_orig1[ens_orig1 %in%  normtissue$Gene ]
normtissue

