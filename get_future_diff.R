

### analyse future diff 

# 1. diff vars
# 2. CHECK IF variables are normal 
# 3. before showing correlations
y=all_diff_variables_prog[1]
sm_pd<-MOFAobjectPD@samples_metadata
sm_pd_diff<-sm_pd[,c(all_diff_variables_prog, 'PATNO_EVENT_ID')]

sm_pd_diff<-apply(sm_pd_diff,2,as.numeric)

shapiro_diff<-apply(sm_pd_diff,2,function(x){
                          if (!all(is.na(x)))
                        shapiro.test(x)$p.value>0.05
                              }
                    )
shapiro_diff$NP2PTOT_diff_V13_V14
shapiro_diff$NP3TOT_diff_V13_V14 ### this one is normal shapiro is true




mydf <- reshape::melt(sm_pd_diff, id.vars=c('PATNO_EVENT_ID'))



library('ggplot2')


p<-ggplot(data =mydf, aes_string(x='value') )+
  geom_density()+ 
  facet_wrap(~variable, scales='free')
  
p
graphics.off()
