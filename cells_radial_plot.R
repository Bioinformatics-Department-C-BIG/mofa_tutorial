
library(reshape2)

cells<-read.csv2('../radar_reshaped.txt',sep = '\t')
cells<-cells[order(cells$Betweenness.C., decreasing = TRUE),]
cells$var.order<-order(cells$Betweenness.C., decreasing = TRUE)


cells_r<-melt(cells,id.vars = c('Cells', 'var.order'))
cells_r$value<-as.numeric(cells_r$value)


# Create a label for the plot in 2 lines to fit by adding \n
cells_r$labels<-paste0(sapply(strsplit(as.character(cells_r$Cells), " "), paste, collapse='\n'))


library(ggplot2)
ggplot(cells_r, aes(y =value , x = reorder(Cells, var.order), 
                group = variable, colour = variable)) + coord_polar(start=0) +
  geom_point() + geom_path(size=1) + 
  theme_bw(base_line_size = 1.2) +
  theme(text=element_text(size=16))+
  theme(legend.position = 'bottom')+
  scale_y_continuous(breaks = seq(0, 4, by = 1), position = 'left') + 
  scale_x_discrete(labels=cells_r$labels)+
  labs(x = NULL, y=NULL)

ggsave('polar.png', height=5.5)



