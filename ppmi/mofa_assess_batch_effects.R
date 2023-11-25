



######### MOFA object covariates #############






covariates=c('SITE', 'Plate', 'NHY', 'NP2PTOT_LOG', 'AGE_SCALED', 'SEX', 
             'Usable.Bases....', 'updrs2_score',  'scopa')

correlate_factors_with_covariates_categ(MOFAobjectPD, covariates=covariates)



kruskal.test(sm_pd$SITE,sm_pd$Plate)
########## boxplots by site and mds2 updrs ####

hist(sm_pd$NP2PTOT_LOG)
sm_pd$SITE=as.factor(sm_pd$SITE)
sm_pd$Plate=as.factor(sm_pd$Plate)

x='Plate'

ggplot(sm_pd, aes_string(y='NP2PTOT',x=x , fill='SITE'))+
  geom_boxplot(aes_string(y='NP2PTOT', fill='SITE'))


fs=c(2,7); color_by='Plate'
MOFAobject@samples_metadata$Plate<-as.numeric(MOFAobject@samples_metadata$Plate)
pf<-plot_factors(MOFAobject, 
                 factors = fs, 
                 #shape_by=color_by,
                 color_by = color_by,
                 show_missing = FALSE 
)


pf

