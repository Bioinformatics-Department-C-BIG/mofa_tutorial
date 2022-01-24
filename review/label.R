all<-read.csv2('review/overlap.txt')
SC<-read.csv2('review/SC.txt')
Method<-read.csv2('review/Method.txt')
Reviews<-read.csv2('review/Reviews.txt')

all[all[,1] %in% Method[,1],'Type']<-'Method'
all[all[,1] %in% Reviews[,1],'Type']<-'Review'

write.csv2(all, 'all_labels.txt')

# overlap between new list

mylist<-read.csv2('review/mylist_601.txt')
new_and_anaysis<-read.csv2('review/metananalysis.txt')


mylist[mylist[,1] %in% new_and_anaysis[,1],'Type']<-'meta-analysis'

write.csv2(mylist[,2], 'review/in_metaanalysis.txt')

