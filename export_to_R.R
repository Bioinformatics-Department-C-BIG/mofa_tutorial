

  tool_obj[tool_obj$method_orig=='pls',]

all_tools<-read_excel('E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/all tools.xlsx')


all_tools_new<-as.data.frame(merge(all_tools[c('method_orig')], tool_obj[c('method_orig', 'objectives_concat')], by='method_orig'))

df1 = all_tools[c('method_orig')]
df2=tool_obj[c('method_orig', 'Var1',  'objectives_concat')]
all_tools_new= df1 %>% left_join(df2,by="method_orig")


write.table(all_tools_new, file = "review/output/tools_objectives_concatenated.txt", sep='\t', row.names = FALSE, quote = FALSE)

# Only output the ones that I described

all_tools_to_write<-all_tools[!is.na(all_tools$Refined),]
all_tools_to_write<-all_tools_to_write[
  with(all_tools_to_write, order(Var1, objectives_concat, method_orig,decreasing = TRUE)),
]


all_tools_to_write$x<-'\\\\'

write.table(all_tools_to_write[c('name display', 'Var1','objectives_concat','x' )], file = "review/output/tools_objectives_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

data<-all_tools_to_write[c('name display', 'Var1','objectives_concat')]




library(RColorBrewer)


library('flextable')


data$Var1=as.factor(data$Var1)

ft1=flextable(data) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))%>%
  highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("Var1"))

ft1

save_as_html( 'Tools'= ft1, 
              path = 'plots/tools_table.html')

save_as_image( ft1, 
               path = 'plots/tools_table.png')





