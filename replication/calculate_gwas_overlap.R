



all_genes_mofa<-sample(rownames(MOFAobject@data$RNA[[1]]), )




#percs=c(0.001, 0.005,0.009,0.01,0.05,0.09,0.1,0.15,0.2)
percs<-c(2:10 %o% 0.1^(2:3))

percs
nfound=c()
nfound_random=c()


fact_coh<-which(cors[,'COHORT']>0)
fact_coh
score='NP3TOT'
score = NULL

score='NP3TOT_LOG'

score='NP2PTOT_LOG'

score = 'COHORT'

score='NP3TOT_LOG'

fact=get_factors_for_metric(score)
fact
#fact=fact[!fact%in% c(13) ]

if (is.null(score)){
    fact = 1:34
} else if (score=='COHORT'){
    fact = fact_coh
}


#fact = intersect(fact, fact_coh)
#fact
for (top_fr in percs ){
    print(paste('top_fr',top_fr))
    top_rnas<-concatenate_top_features(MOFAobject, factors = fact,  top_fr=top_fr, view='RNA')
    top_rnas

    unique(top_rnas$feature)
    top_mofa_weights<-unique(get_symbols_vector(top_rnas$feature))

    number_perc<-length(top_mofa_weights)


    gwas<-read.table('ppmi/gwas_genes_kim2024.txt')

    length(gwas$V1)
    any(top_mofa_weights %in% gwas$V1)
    which(top_mofa_weights %in% gwas$V1)
    print(length(which(top_mofa_weights %in% gwas$V1)))
    nfound_perc=length(which(top_mofa_weights %in% gwas$V1))


    ### assess the random 

    all_random = sapply(1:30,function(i){
        random_genes<-unique(get_symbols_vector(sample(all_genes_mofa,number_perc)))

        nfound_random_perc = length(which(random_genes %in% gwas$V1))
        return(nfound_random_perc)
    })

    nfound_random_perc_average<-mean(all_random)

    nfound=c(nfound,nfound_perc )
    nfound_random=c(nfound_random,nfound_random_perc_average )

}




nfound_random
fact_s=paste0(fact, collapse='_')
fact_s

suppressWarnings(dir.create(paste0(outdir,'/gwas_overlap/')))
graphics.off()
jpeg(paste0(outdir,'/gwas_overlap/',score,fact_s,  '.jpeg'), res=100)
plot(percs, nfound, col='blue')
points(percs, nfound_random, col='red')
dev.off()


plot(percs, nfound, col='blue')
points(percs, nfound_random, col='red')





