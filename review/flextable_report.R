
library(RColorBrewer)





data$Var1=as.factor(data$Var1)

ft1=flextable(data) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))%>%
  highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("Var1"))



save_as_html( 'Tools'= ft1, 
  path = 'plots/tools_table.html')

save_as_image( ft1, 
              path = 'plots/tools_table.png')

