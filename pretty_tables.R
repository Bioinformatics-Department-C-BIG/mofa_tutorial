
library(tidyverse)
library(glue)
library(gt)



# 
# install.packages('reactable')
library(htmltools)
# 
# install.packages('reactablefmtr')




table_data<-data %>% reactable(columns = list(Var1 = colDef( cell =
                                                               color_tiles(., number_fmt = scales::col_factor(Var1),
                                                                           colors = viridis::viridis(5))
)))

my_color_pal = c('#e5f5e0', '#a1d99b', '#31a354')

reactable(
  data, 
  columns=list( 
    Var1 = colDef(
      style = color_scales(data))
  )
)


save_reactable_test(table_data, "table.png")


library(viridis)




flextable(data) %>%
  bg(bg = scales::col_numeric(palette = "viridis", domain = c(0, 1)))

