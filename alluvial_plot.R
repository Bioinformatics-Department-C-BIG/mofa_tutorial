
################ Section 5: ALLUVIAL 
# This depends on inheriting the data frame 'new' from literature_review.R
# Run it with the appropriate x_group before this script. 

#install.packages('alluvial')
#install.packages('ggalluvial')
#install.packages('ggsankey')
library('ggalluvial')
library('alluvial')

#library('ggsankey')
library('dplyr')
library('tidyr')

# Switch this to print both
cancer_filter=c("no")

new2<-new %>% 
  mutate(Data=strsplit(Data, ',|\r|\n' ) )%>%
  unnest(Data) 
new2$objective
# for data-disease run literature_review, for objective run combined
axis1='Data'
axis2='objective'

if (axis2 == 'disease_group'){
  new2$disease_group<-group_disease(new2, 'Disease')$Disease
  new2<-new2[new2$disease_group %in% df_most_common_disease$disease_group,]}
new2<-new2 %>% filter( (!!sym(axis2)) %in% most_common_groups)
new2<-new2 %>% filter( ! (!!sym(axis2)) %in% remove_objectives)
new2<-group_objectives(new2, 'objective')

new2<-new2[!is.na(new2[axis2]),]
new2<-new2[!is.na(new2[axis1]),]

new2<-new2 %>% filter(Cancer %in% cancer_filter)


new2$Data<-tolower(trimws(new2$Data))
new2<-new2 %>% filter(Data %in% tolower(level2))

levels(as.factor(new2$Data))



axis2='objective'

counts<-new2 %>%
  count(objective, method)
new2
n_cutoff<-2

df<-counts
counts<-new2 %>% 
  count(objective, method)




##### Count most common objective
unique_ids<-new2 %>% 
  group_by_at('objective') %>%
  count('PMID')
NROW(unique_ids)



# Count unique objective-study pairs 
new2 %>% 
  group_by_at(c('objective', 'PMID')) %>%
  count(c('objective', 'PMID'))%>%
  group_by_at(c('objective')) %>%
  count('PMID')




# New implementation with ggalluvial
axis1='Data'
axis2="objective"



counts1 <- new2 %>% 
  group_by_at(c('Data', 'objective')) %>%
  count(c('Data', 'objective')) %>%
  mutate(
     col = objective
  )%>%
  as.data.frame()
n_cutoff=5

#Choose top omics only 
omics_freqs= as.data.frame(aggregate(counts1$freq, list(counts1$Data), FUN=sum) )
top_omics<-omics_freqs[order(omics_freqs$x, decreasing = TRUE),1][1:4]
counts1=counts1[counts1$Data %in% top_omics,]
counts1

# Rename the objectives 

counts1$key_names<-counts1$objective
counts1<-relabel_objectives_short(counts1)
counts1$objective<- counts1$key_names
counts1$labels<-NULL
counts1$col=counts1$objective
  

counts1<-counts1[!(counts1$objective %in% remove_objectives),]


# Updated to filter groups as before
counts1<-counts1%>% 
  group_by('Data') %>%
  #filter(n>n_cutoff)

counts1
 counts<-counts1 %>% ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))


 df<-counts
g<-ggplot(data = df, aes(x = type, stratum = stratum, alluvium = alluvium, y = freq)) +
  geom_flow(aes(fill = col), width = 1/6, color = "darkgray",
            curve_type = "cubic") +
  geom_stratum(color = "grey", width = 1/6) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme(
    panel.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  ) +
  scale_fill_viridis_d()


if (cancer_filter=='no'){
   g<-g+ggtitle(paste0("Other Diseases"))
  }else{
  g<-g+ggtitle(paste0("Cancer"))
  }
g<-g+    theme(plot.title = element_text(hjust = 0.5))

g

fname=paste0('plots/ggalluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png')
ggsave(fname, width = 7, height=7)

