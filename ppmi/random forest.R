


suppressPackageStartupMessages(library(randomForest))
# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
y='NHY'
df$res
MOFAobject@samples_metadata[y ]
df$res
MOFAobject@samples_metadata[y ][[1]]
df$res <- as.factor(as.numeric(MOFAobject@samples_metadata[y ][[1]]  ))
model.ighv <- randomForest(res ~ ., data=df[!is.na(df$res),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)
# Train the model for Trisomy12
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
model.trisomy12 <- randomForest(trisomy12 ~ ., data=df[!is.na(df$trisomy12),], ntree=10)
df$trisomy12 <- NULL

MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)
Plot predictions for IGHV

MOFAobject@samples_metadata$IGHV.pred_logical <- c("True","Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV))+1]

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
