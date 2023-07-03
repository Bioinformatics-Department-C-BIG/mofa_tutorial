
library(ggplot2)

source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/utils.R'))
## run mofa / or just load model 
source(paste0(script_dir,'ppmi/mofa_application_ppmi_all_visits.R'))
library(R.filesets)
library(stats)

library(maxstat )
install.packages('survminer')

library(survminer )



### Variable to predict: 
# 1. where there is change in HY or not? 
# 2. 

### TODO: plot clinical change by factor 
### Detect baseline clinical change by factor.


MOFAobject@samples_metadata$Outcome
se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)

if (filter_common){
  se_filt_BL_common<-se_filt_BL[,se_filt_BL$PATNO %in% common]
  se_filt_V08_common<-se_filt_V08[,se_filt_V08$PATNO %in% common]
}


dim(se_filt_BL_common)
dim(se_filt_V08_common)
### very important: reorder: 
se_filt_BL_common=se_filt_BL_common[,match(se_filt_V08_common$PATNO,se_filt_BL_common$PATNO) ]


change<-se_filt_V08_common$NP3TOT-se_filt_BL_common$NP3TOT; hist(change)
quantile(change, 0.9, na.rm=TRUE)
change_discrete<-ifelse(change>quantile(change, 0.9, na.rm=TRUE), 1, 0); change_discrete
hist(se_filt_V08_common$NP2PTOT-se_filt_BL_common$NP2PTOT)

se_filt_V08_common$change_discrete<-change_discrete
se_filt_V08_common$change<-change
col

colData(se_filt_V08_common)[, c( 'PATNO_EVENT_ID', 'change', 'change_discrete')]
se_filt_V08_common_mfa<-se_filt_V08_common[, match(MOFAobject@samples_metadata$PATNO,se_filt_V08_common$PATNO )]

Z <- get_factors(MOFAobject)[[1]]

ordered_vars<-colData(se_filt_V08_common_mfa)[c('change', 'change_discrete')]
rownames(ordered_vars)<-MOFAobject@samples_metadata$PATNO_EVENT_ID
ordered_vars[15,]
MOFAobject@samples_metadata$PATNO[15]
MOFAobject@samples_metadata<-cbind(MOFAobject@samples_metadata, ordered_vars[,c('change', 'change_discrete')])
Z1=Z[,sel_factors[4]]
MOFAobject@samples_metadata$change


df<-data.frame(cbind(MOFAobject@samples_metadata, Z1=Z1))



hist(Z1)

ggplot(df, aes(x=Z1))+
  geom_density()


ggplot(df, aes(x=Z1, y=change))+
  geom_point()


####


df <- data.frame(
  time = SurvObject[,1], 
  event = change, Z1 = Z[,1]
)





df<-data.frame(cbind(MOFAobject@samples_metadata, Z1=Z[,sel_factors[1]]))
df$Z1






cut <-surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint
fit <- survfit(Surv(time, event) ~ FactorCluster, df)

ggsurvplot(fit, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100,
           legend = "top", legend.labs = c(paste("low LF 1"), paste("high LF 1")),
           xlab = "Time to treatment", ylab="Survival probability (%)", title= "Factor 1"
)$plot