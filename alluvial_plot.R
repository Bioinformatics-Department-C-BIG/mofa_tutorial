
################ Section 5: ALLUVIAL 


#install.packages('alluvial')
#install.packages('ggalluvial')
#install.packages('ggsankey')
library('ggalluvial')
library('alluvial')

library('ggsankey')

cancer_filter=c("yes")
new2<-new %>% 
  mutate(Data=strsplit(Data, ',|\r|\n' ) )%>%
  unnest(Data) 


#new2<-new2 %>% filter(Cancer %in% cancer_filter)


new2$Data<-tolower(trimws(new2$Data))
new2<-new2 %>% filter(Data %in% tolower(level1))

levels(as.factor(new2$Data))


new2<-new2[!is.na(new2$method),]
new2<-new2[!is.na(new2$Data),]
new2<-new2[!is.na(new2$objective),]

axis1='Data'
axis2='objective'

counts<-new2 %>% count(objective, method)

if (cancer_filter == 'yes'){n_cutoff=3} else {n_cutoff=4}
counts<-counts%>% filter(n>n_cutoff)


df<-counts
counts<-new2 %>% count(objective, method)
counts<-counts%>% filter(n>n_cutoff)




# New implementation with ggalluvial
axis1='Data'
axis2="objective"
counts <- new2 %>% 
  count( Data, objective) %>% 
  mutate(
    col = objective
  ) %>%
  ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))

df<-counts %>% filter(n>n_cutoff)
# df<-counts
ggplot(data = df, aes(x = type, stratum = stratum, alluvium = alluvium, y = n)) +
  # geom_lode(width = 1/6) +
  geom_flow(aes(fill = col), width = 1/6, color = "darkgray",
            curve_type = "cubic") +
  # geom_alluvium(aes(fill = stratum)) +
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
  scale_fill_viridis_d()+
  ggtitle(paste0("Multi omics objectives, Cancer = ", cancer_filter))


ggsave(paste0('plots/ggalluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png'), width = 7, height=6)
