

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

all_genes_mofa<-sample(rownames(MOFAobject@data$RNA[[1]]), )
se_rnas
all_genes_mofa<-sample(rownames(MOFAobject@data$RNA[[1]]), )

all_genes_mofa = rownames(se_rnas)

se_rnas
all_genes_mofa

#percs=c(0.001, 0.005,0.009,0.01,0.05,0.09,0.1,0.15,0.2)
percs<-c(2:10 %o% 0.1^(2:3))

percs<-c(0.05)
percs<-c(0.2)

nfound=c()
nfound_random=c()


fact_coh<-which(cors[,'COHORT']>0)
fact_coh
score='NP3TOT'
score = NULL

score='NP3TOT_LOG'




score='NP2PTOT_LOG'
score='NP3TOT_LOG'
score='sft'

fact=get_factors_for_metric(score)
fact
#score = 'COHORT'

#fact=fact[!fact%in% c(13) ]

if (is.null(score)){
    fact = 1:34
} else if (score=='COHORT'){
    fact = fact_coh
}


#fact = intersect(fact, fact_coh)
#fact
use_jaccard = TRUE
view='RNA'
for (top_fr in percs ){
    print(paste('top_fr',top_fr))
    top_rnas<-concatenate_top_features(MOFAobject, factors = fact,  top_fr=top_fr, view=view)
    top_rnas

    unique(top_rnas$feature)
    top_mofa_weights<-unique(convert_to_gene_symbol(top_rnas$feature, view=view))

    number_perc<-length(top_mofa_weights)
    print(top_mofa_weights)


    gwas<-read.table('ppmi/gwas_genes_kim2024.txt')

     gwas<-read.table('ppmi/gwas_genes_kim2024.txt');     known_genes<-gwas$V1
   disgenet<-read.csv2('/data8TB/efiath/git/mofa_tutorial/replication/C0030567_disease_gda_summary-PD.txt', sep='\t'); 
   disgenet<-disgenet[disgenet$Score_gda>0.1,]
   disgenet

    known_genes= disgenet$Gene
       known_genes= gwas$V1
    print(length(known_genes))


   

    length(gwas$V1)
    any(top_mofa_weights %in%  known_genes)
    which(top_mofa_weights %in%  known_genes)
    print(length(which(top_mofa_weights %in%  known_genes)))
    nfound_perc=length(which(top_mofa_weights %in%  known_genes))
    
    jo_index = jaccard( known_genes,top_mofa_weights)
    jo_index

    inter<-intersect(top_mofa_weights, known_genes)
    print(top_rnas[top_rnas$feature %in% inter,])
    if (use_jaccard){
        nfound_perc = jo_index

    }


    ### assess the random 

    all_random = sapply(1:5,function(i){
        random_genes<-unique(get_symbols_vector(sample(all_genes_mofa,number_perc)))

        nfound_random_perc = length(which(random_genes %in% known_genes))
        print(length(random_genes))
        print(length(known_genes))
        jo_index_random = jaccard( known_genes,random_genes)

        if (use_jaccard){
            nfound_random_perc = jo_index_random
        }




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
print(paste0(outdir,'/gwas_overlap/',score,fact_s,'_',use_jaccard,  '.jpeg'))
jpeg(paste0(outdir,'/gwas_overlap/',score,fact_s,'_',use_jaccard,  '.jpeg'), res=100)
plot(percs, nfound, col='blue')
points(percs, nfound_random, col='red')
dev.off()


plot(percs, nfound, col='blue')
points(percs, nfound_random, col='red')

#print(paste(score, round(nfound[5], digits=3), round(nfound_random[5], digits=3)))
print(paste(score, round(nfound, digits=3), round(nfound_random, digits=3)))






