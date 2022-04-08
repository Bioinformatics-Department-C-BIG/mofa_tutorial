


#### 

  tool_obj[tool_obj$method_orig=='pls',]

all_tools<-read_excel(paste0(data_int_dir,'/all tools.xlsx'))


all_tools_new<-as.data.frame(merge(all_tools[c('method_orig')], tool_obj[c('method_orig', 'objectives_concat')], by='method_orig'))

df1 = all_tools
df2=tool_obj[c('method_orig', 'Var1',  'objectives_concat')]
all_tools_new= df1 %>% full_join(df2,by="method_orig")

all_tools_new2<-sapply(all_tools_new,as.character)
all_tools_new2[is.na(all_tools_new2)]<-""


write.table(all_tools_new2, file = "review/output/tools_objectives_concatenated.txt", sep='\t', row.names = FALSE, quote = FALSE)

# Only output the ones that I described

all_tools_to_write<-all_tools[!is.na(all_tools$Refined),]
all_tools_to_write<-all_tools_to_write[
  with(all_tools_to_write, order(Var1, objectives_concat, method_orig,decreasing = TRUE)),
]


all_tools_to_write$x<-'\\\\'


# todo: add column to all_tools ? 
write.table(all_tools_to_write[c('name display', 'Var1','objectives_concat','Tool/Method' ,'x' )], file = "review/output/tools_objectives_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

data<-all_tools_to_write[c('name display', 'Var1','objectives_concat', 'Tool/Method')]


# 
# 
# footnote(
#   x,
#   i = NULL,
#   j = NULL,
#   value,
#   ref_symbols = NULL,
#   part = "body",
#   inline = FALSE,
#   sep = "; "
# )



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



data2<-obj_tool_wide[,c('objective', 'Var1', 'no', 'yes')]
data2<-obj_tool_wide[,c('objective', 'Var1', 'no', 'yes')]

ft2=flextable(data2) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))%>%
  highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("no", "yes"), source=c("Var1"))


ft2<-void(ft2, j='Var1', part='body')
ft2<-merge_v(ft2, j = 'objective', target = c('objective'), part = "body", combine = FALSE)

ft2
save_as_image( ft2, 
               path = 'plots/objective_tools_table.png')


