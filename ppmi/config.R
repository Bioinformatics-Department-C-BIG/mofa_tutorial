
MIN_COUNT_G=100
MIN_COUNT_M=10
VISIT='BL'

VISIT=c('V04')
VISIT=('BL')
VISIT=('BL')




TOP_GN=0.1
TOP_MN=0.5


sel_coh=c(1,4)


sel_coh <- c(1,2)


sel_coh=c(1,2);





VISIT=c( 'V08')
VISIT=c('V06')

VISIT_S=paste(VISIT,sep='_',collapse='-')
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')

g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_')
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  )

