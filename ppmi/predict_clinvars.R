
library(GGally)

y<-MOFAobject@samples_metadata$NP3_TOT
MOFAobject@samples_metadata$SEX

w_fs<-get_weights(MOFAobject, view='RNA', factors = c(4))[[1]]
selected_feats<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.995),])
x_rnas<-MOFAobject@data$RNA[[1]][selected_feats, ]

w_fs<-get_weights(MOFAobject, view='miRNA', factors = c(3))[[1]]
selected_feats2<-names(w_fs[abs(w_fs)>quantile(abs(w_fs), 0.995),])


x_prots<-MOFAobject@data$miRNA[[1]][selected_feats2, ]

x_prots
x_feats<-rbind(x_rnas, x_prots)
length(selected_feats)
####


### SELECT FACTORS ####
Z <- get_factors(MOFAobject)[[1]]
x_factors<-data.frame(cbind(Z[,c(sel_factors, 15)], y))
x_factors
GGally::ggpairs(x_factors)
x_vars<-cbind(Z[,c(sel_factors, 15)],MOFAobject@samples_metadata[, c('AGE', 'SEX')]   )
x_vars<-cbind(Z[,c(sel_factors, 15)])
x_vars<-cbind(Z[,c(sel_factors, 15)],MOFAobject@samples_metadata[, c('AGE')]   )
x_vars<-cbind(Z[,c(sel_factors, 15)]  )


#### SELECT MOLECULES ####
x_vars<-cbind(t(MOFAobject@data$RNA[[1]][selected_feats, ]))
x_vars<-cbind(t(x_feats), MOFAobject@samples_metadata[, c('AGE')] )

x_vars<-cbind(t(x_feats), MOFAobject@samples_metadata[, c('AGE')] )


dim(x_vars)
dim(x_vars)

y
fit <- lm(y ~ Z)
fit <- lm(y ~ x_vars)

fit$coefficients
fit_an<-anova(fit)
fit_an
f <- summary(fit)
f$r.squared
f

sel_factors
