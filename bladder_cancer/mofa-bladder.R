if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")
devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
browseVignettes("MOFA2")
BiocManager::install("MOFAdata")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

row.names(X1_raw)
data("CLL_data")


data = list(RNA=as.X1_raw, 
            proteomics = X2_raw )
MOFAobject <- create_mofa(data)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

model_opts$likelihoods


#### Training 
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts


##### Prepare MOFA 
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

##### slots 
MOFAobject <- runMOFA(MOFAobject)

MOFAobject <- readRDS(url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds"))

slotNames(MOFAobject)
names(MOFAobject@data)

dim(MOFAobject@data$Drugs$group1)

names(MOFAobject@expectations)

dim(MOFAobject@expectations$Z$group1)

# mrna weight matrix 
dim(MOFAobject@expectations$W$mRNA)


samples_metadata(MOFAobject) <- CLL_metadata
plot_factor_cor(MOFAobject)


plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=15)





plot_variance_explained(MOFAobject, plot_total = T)[[2]]


#### Association 
### Test the association between MOFA factors and gender, survival outcome 
### and age 

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)



#### Factor analysis 
### like principal component analysis 

  plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)



### Plot feature weights 
### Feature weights for mutations 

plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


###Plot gene weights for MRNA expression #

plot_weights(MOFAobject, 
             view = "mRNA", 
             factor = 1, 
             nfeatures = 10
)


utils::data(reactomeGS)


plot_weights(MOFAobject, 
             view = "Mutations", 
             factor = 3, 
             nfeatures = 10,
             abs = F
)


##### Inspection of combinations of Factors

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV",
                  shape_by = "trisomy12",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)



###### Heatmaps 

  plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 5,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)



##### Predicton of cilincal subgroups 
install.packages('randomForest')
suppressPackageStartupMessages(library(randomForest))


df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(IGHV ~ ., data=df[!is.na(df$IGHV),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)



# Train the model for Trisomy12
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
model.trisomy12 <- randomForest(trisomy12 ~ ., data=df[!is.na(df$trisomy12),], ntree=10)
df$trisomy12 <- NULL

MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)


OFAobject@samples_metadata$IGHV.pred_logical <- c("True","Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV))+1]

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV.pred",
                  shape_by = "IGHV.pred_logical",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)


###### GSEA 
library(MOFAdata)

utils::data(reactomeGS)

head((reactomeGS))

head(colnames(reactomeGS))

# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)



# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "negative"
)


MOFAobject@data$Mutations
res.negative

theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

names(res.positive)

plot_enrichment_heatmap(res.positive)
ggsave(paste0('Enrichment_heatmap_positive','.png'), width = 9, height=4, dpi=100)


plot_enrichment_heatmap(res.negative)
ggsave(paste0('Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=100)


factor_to_plot=9
plot_enrichment(res.positive, factor = factor_to_plot, max.pathways = 15)
ggsave(paste0('GSEA_factor_',factor_to_plot,'.png'), width = 9, height=4, dpi=100)

factor_to_plot=9
plot_enrichment(res.negative, factor = factor_to_plot, max.pathways = 15)
ggsave(paste0('GSEA_factor_neg_',factor_to_plot,'.png'), width = 9, height=4, dpi=100)


##### ASSOCIATE FACTORS WITH SURVIVAL 
library(survival)
library(survminer)


SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z) 
fit


##### Plot hazard ratios 
s <- summary(fit)
coef <- s[["coefficients"]]

df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p      = coef[,"Pr(>|z|)"], 
  coef   = coef[,"exp(coef)"], 
  lower  = s[["conf.int"]][,"lower .95"], 
  higher = s[["conf.int"]][,"upper .95"]
)

ggplot(df, aes(x=factor, y=coef, ymin=lower, ymax=higher)) +
  geom_pointrange( col='#619CFF') + 
  coord_flip() +
  scale_x_discrete() + 
  labs(y="Hazard Ratio", x="") + 
  geom_hline(aes(yintercept=1), linetype="dotted") +
  theme_bw()

