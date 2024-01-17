
#########
# TODO: make a run with all times and plot like in molecular trajectories 





cors_all[, 'COHORT']

f=3
fgs<-as.data.frame(get_factors(MOFAobject, factors=f)[[1]])
sm<-samples_metadata(object = MOFAobject)
#rownames(fgs1)<-

### Match rownames to patient numbers 
y_cs='NP2PTOT'
fgs$EVENT_ID<-sm$EVENT_ID
fgs$NP3_TOT<-sm$NP3_TOT
fgs[,y_cs ]<-sm[, y_cs]

fgs$PATNO<-sm$PATNO
fgs$COHORT<-sm$COHORT


fgs


intersect(rownames(fgs$BL), rownames(fgs$V08))


mf_time_m<-reshape2::melt(fgs)
MOFAobject@samples_metadata[[1]]



y=paste0('Factor', f)
ggplot(fgs)+ 
  geom_point( aes_string(x='EVENT_ID', fill='PATNO', y=y) ) + 
  geom_line(aes_string( y=y, x='EVENT_ID', group='PATNO', color='PATNO'))+
  geom_smooth( aes_string(y=y,  x='EVENT_ID'))+
  facet_wrap(~COHORT, nrow=2)


ggsave(paste0(outdir, '/trajectories/factor', f, '.png'))




#### Ratios ####

length(fgs[fgs$EVENT_ID=='V08',]$PATNO)
length(fgs[fgs$EVENT_ID=='BL',]$PATNO)

fgs_BL_new<-fgs[fgs$EVENT_ID=='BL',]
fgs_V08_new<-fgs[fgs$EVENT_ID=='V08',]


fgs_BL_new<-fgs[fgs$EVENT_ID=='BL',]


common<-intersect(fgs_BL_new$PATNO, fgs_V08_new$PATNO)
fgs_BL_new=fgs_BL_new[fgs_BL_new$PATNO%in%common,]
fgs_V08_new=fgs_V08_new[fgs_V08_new$PATNO%in%common,]

fgs_V08_new=fgs_V08_new[match(fgs_BL_new$PATNO, fgs_V08_new$PATNO),]
fgs_BL_new$PATNO

fgs_V08_new[, y]/fgs_BL_new[, y]
v08_f<-fgs_V08_new[, y]
bl_f<-fgs_BL_new[, y]


min_all<-min(c(v08_f,bl_f ))
max_all<-max(c(v08_f,bl_f ))
scale_values <- function(x, min_x, max_x){
                    (x-min_x)/(max_x-min_x)
                    }
hist(scale_values(v08_f, min_x=min_all, max_x=max_all))

ratios<-log2(scale_values(v08_f, min_x=min_all, max_x=max_all) / scale_values(bl_f, min_x=min_all, max_x=max_all))
ratios<-ratios

outl<-ratios< (-4)
ratios_np3<- log2(fgs_V08_new[, y_cs] /fgs_BL_new[, y_cs])

hist(ratios)
hist(ratios_np3)

plot(ratios, ratios_np3 )
rcs<-data.frame(ratios_np3=ratios_np3, ratios=ratios)

ggplot(rcs[!outl,], aes(x=ratios, y=ratios_np3))+
  geom_smooth()+geom_point()+
  labs(x='Factor log2FC')


cor(ratios)

cors_all[, 'COHORT']


fgs<-get_factors(MOFAobject, factors=7)
sm<-samples_metadata(object = MOFAobject)
smbl<-sm[sm$group=='BL',]
#rownames(fgs1)<-

### Match rownames to patient numbers 
rownames(fgs$BL)<-smbl$PATNO
rownames(fgs$V04)<-sm[sm$group=='V04',]$PATNO
rownames(fgs$V06)<-sm[sm$group=='V06',]$PATNO
rownames(fgs$V08)<-sm[sm$group=='V08',]$PATNO


intersect(rownames(fgs$BL), rownames(fgs$V08))

mf_time<-merge(fgs$BL,fgs$V08, by='row.names')
mf_time_m<-reshape2::melt(mf_time)
MOFAobject@samples_metadata[[1]]



ggplot(mf_time_m)+ 
  geom_point( aes(x=variable, fill=Row.names, y=value) ) + 
  geom_line(aes( y=value, x=variable, group=Row.names))

