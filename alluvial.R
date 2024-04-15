#install.packages('ggsankey')
library(nycflights13)
library('ggsankey')
library('ggplot2')
#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb) 



top_dest <- flights %>% 
  count(dest) %>% 
  top_n(5, n) %>% 
  pull(dest)

top_carrier <- flights %>% 
  filter(dest %in% top_dest) %>% 
  count(carrier) %>% 
  top_n(4, n) %>% 
  pull(carrier)

fly <- flights %>% 
  filter(dest %in% top_dest & carrier %in% top_carrier)


fly <- flights %>% 
  filter(dest %in% top_dest & carrier %in% top_carrier) %>%
  ggsankey::make_long(origin, carrier, dest)




fly <- flights %>% 
  filter(dest %in% top_dest & carrier %in% top_carrier) %>%
  count(origin, carrier, dest) %>% 
  mutate(
    origin = fct_relevel(as.factor(origin), c("LGA", "EWR","JFK")),
    col = origin
  ) %>%
  ggalluvial::to_lodes_form(key = type, axes = c("origin", "carrier", "dest"))
ggplot(data = fly, aes(x = type, stratum = stratum, alluvium = alluvium, y = n)) +
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
  scale_fill_viridis_d()

