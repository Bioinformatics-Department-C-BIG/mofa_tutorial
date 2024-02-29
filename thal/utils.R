

preprocess_clinical_vars<-function(MOFAobject){
    sm<-MOFAobject@samples_metadata

    total_vols<-colnames(sm)[grep('Total|Number', colnames(sm))]
    sm[,total_vols]<-sapply(sm[,total_vols], as.numeric)
    sm[,total_vols]
    sm$STAGE<-as.factor(sm$STAGE)
    return(sm)
}


