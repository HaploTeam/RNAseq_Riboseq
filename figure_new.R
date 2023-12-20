####_####
####Main figures ####
####_####

#### Figure 1 ####
#####A#####
x= cor(mean_norm_Ribo, method = 's')
x=melt(x)
temp=x$value
x= cor(mean_norm_RNA, method = 's')
x=melt(x)
temp= median(c(temp,x$value))
x= cor(mean_norm_Ribo, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(temp),lim =c(0.4,1),name='Correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)+
  theme_minimal()+
  xlab('')+
  ylab('')

x= cor(mean_norm_RNA, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(temp),lim =c(0.4,1), name='Correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)+
  theme_minimal()+
  xlab('')+
  ylab('')

#####B#####
x= cor(mean_norm_RNA, method = 's')[lower.tri(cor(mean_norm_RNA, method = 's'))]
temp= cor(mean_norm_Ribo, method = 's')[lower.tri(cor(mean_norm_Ribo, method = 's'))]
temp = data.frame(Ribo= temp, RNA=x)
ggplot(temp, aes(RNA, Ribo))+
  geom_point()+
  theme_classic()+
  stat_cor()+
  xlab('RNA correlation indexes')+
  ylab('Ribo correlation indexes')

#####C#####

DE_gene = unique(one_vs_all_Ribo[one_vs_all_Ribo$is_DE,"gene"])
goterms <- Term(GOTERM)
x= goterms
x= goterms[grep('metabolic',x )]
x= names(x)
x= lapply(geneID2GO,function(i){
  i<<-i
  length(intersect(i,x))>0
})
x= unlist(x)
metabo_gene = names(geneID2GO)[x]

a = sum(unique(DE_gene)%in%metabo_gene)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_Ribo)%in%metabo_gene)-a
d = length(rownames(norm_Ribo))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)$estimate
fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate
x= data.frame(Odd_Ratio=c(fisher.test(e)$estimate,
                          fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate)
)
x$type= c('Metabolism','Other')
p1= ggplot(x, aes(type,-1+Odd_Ratio, fill = type))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1,linetype='dashed',col='grey', size=1.5)+
  geom_vline(xintercept = 2,linetype='dashed',col='grey',size=1.5)+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('gray23','gray71'))+
  coord_flip()+
  theme_minimal()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5, 1), labels = c(0,0.5,1,1.5,2))+
  xlab('')+
  ylab('Odds Ratios')+
  theme(legend.position = 'none')

colnames(e)= c('Metabolic', 'Other')
rownames(e)= c('DE', 'Not DE')
e = rbind(e,Total= colSums(e))

e= melt(e)

x=e[e$Var1=='DE',]
x$per= x$value/sum(x$value)*100
a=e[e$Var1=='Not DE',]
a$per= a$value/sum(a$value)*100
b=e[e$Var1=='Total',]
b$per= b$value/sum(b$value)*100

x = rbind(a,x)
x=rbind(x,b)


p2= ggplot(x[!x$Var1=='Total',], aes(fill = Var2, y=per, x= Var1))+
  geom_bar(stat = 'identity')+
  coord_polar(theta = "y", start = 0, direction = 1, clip = "on")+
  scale_fill_manual(values = c('gray23','gray71'))+
  theme_void()+
  geom_bar(data= x[x$Var1=='Total',], aes(fill = Var2, y=per, x= 2.6),stat = 'identity', width=0.2)


grid.arrange(p2,p1,ncol=1)

#####D#####

DE_gene = unique(one_vs_all_RNA[one_vs_all_RNA$is_DE,"gene"])
goterms <- Term(GOTERM)
x= goterms
x= goterms[grep('metabolic',x )]
x= names(x)
x= lapply(geneID2GO,function(i){
  i<<-i
  length(intersect(i,x))>0
})
x= unlist(x)
metabo_gene = names(geneID2GO)[x]

a = sum(unique(DE_gene)%in%metabo_gene)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_RNA)%in%metabo_gene)-a
d = length(rownames(norm_Ribo))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)$estimate
fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate
x= data.frame(Odd_Ratio=c(fisher.test(e)$estimate,
                          fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate)
)
x$type= c('Metabolism','Other')
p1= ggplot(x, aes(type,-1+Odd_Ratio, fill = type))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 1,linetype='dashed',col='grey', size=1.5)+
  geom_vline(xintercept = 2,linetype='dashed',col='grey',size=1.5)+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('gray23','gray71'))+
  coord_flip()+
  theme_minimal()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5, 1), labels = c(0,0.5,1,1.5,2))+
  xlab('')+
  ylab('Odds Ratios')+
  theme(legend.position = 'none')

colnames(e)= c('Metabolic', 'Other')
rownames(e)= c('DE', 'Not DE')
e = rbind(e,Total= colSums(e))

e= melt(e)

x=e[e$Var1=='DE',]
x$per= x$value/sum(x$value)*100
a=e[e$Var1=='Not DE',]
a$per= a$value/sum(a$value)*100
b=e[e$Var1=='Total',]
b$per= b$value/sum(b$value)*100

x = rbind(a,x)
x=rbind(x,b)


p2= ggplot(x[!x$Var1=='Total',], aes(fill = Var2, y=per, x= Var1))+
  geom_bar(stat = 'identity')+
  coord_polar(theta = "y", start = 0, direction = 1, clip = "on")+
  scale_fill_manual(values = c('gray23','gray71'))+
  theme_void()+
  geom_bar(data= x[x$Var1=='Total',], aes(fill = Var2, y=per, x= 2.6),stat = 'identity', width=0.2)


grid.arrange(p2,p1,ncol=1)


#####E#####

x = lapply(unique(one_vs_all_Ribo$isolate),function(i){
  temp= i
  i = one_vs_all_Ribo[one_vs_all_Ribo$isolate==i&one_vs_all_Ribo$is_DE,"gene"]
  
  a = sum(unique(i)%in%unlist(comp_new))
  b = length(unique(i))-a
  c = sum(rownames(norm_Ribo)%in%unlist(comp_new))-a
  d = length(rownames(norm_Ribo))-a-b-c
  e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
  x= fisher.test(e)
  data.frame(OR=x$estimate, p_val=x$p.value,isolate=temp)
  
  
})
x=do.call(rbind ,x)
x$data='Ribo'
x$type= 'Complex'
deg_ess_comp = x
x = lapply(unique(one_vs_all_Ribo$isolate),function(i){
  temp= i
  i = one_vs_all_Ribo[one_vs_all_Ribo$isolate==i&one_vs_all_Ribo$is_DE,"gene"]
  
  a = sum(unique(i)%in%ess)
  b = length(unique(i))-a
  c = sum(rownames(norm_Ribo)%in%ess)-a
  d = length(rownames(norm_Ribo))-a-b-c
  e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
  x= fisher.test(e)
  data.frame(OR=x$estimate, p_val=x$p.value,isolate=temp)
  
  
})
x=do.call(rbind ,x)
x$data='Ribo'
x$type= 'Essential'
deg_ess_comp=rbind(deg_ess_comp,x)

x = lapply(unique(one_vs_all_RNA$isolate),function(i){
  temp= i
  i = one_vs_all_RNA[one_vs_all_RNA$isolate==i&one_vs_all_RNA$is_DE,"gene"]
  
  a = sum(unique(i)%in%unlist(comp_new))
  b = length(unique(i))-a
  c = sum(rownames(norm_RNA)%in%unlist(comp_new))-a
  d = length(rownames(norm_Ribo))-a-b-c
  e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
  x= fisher.test(e)
  data.frame(OR=x$estimate, p_val=x$p.value,isolate=temp)
  
  
})
x=do.call(rbind ,x)
x$data='RNA'
x$type= 'Complex'
deg_ess_comp=rbind(deg_ess_comp,x)

x = lapply(unique(one_vs_all_RNA$isolate),function(i){
  temp= i
  i = one_vs_all_RNA[one_vs_all_RNA$isolate==i&one_vs_all_RNA$is_DE,"gene"]
  
  a = sum(unique(i)%in%ess)
  b = length(unique(i))-a
  c = sum(rownames(norm_RNA)%in%ess)-a
  d = length(rownames(norm_Ribo))-a-b-c
  e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
  x= fisher.test(e)
  data.frame(OR=x$estimate, p_val=x$p.value,isolate=temp)
  
  
})
x=do.call(rbind ,x)
x$data='RNA'
x$type= 'Essential'
deg_ess_comp=rbind(deg_ess_comp,x)
deg_ess_comp$corrected = -1+deg_ess_comp$OR
deg_ess_comp$p_val_corr= p.adjust(deg_ess_comp$p_val,method = 'fdr')
deg_ess_comp$data_plot= deg_ess_comp$data
deg_ess_comp[deg_ess_comp$p_val_corr>0.05,"data_plot"]=paste0(deg_ess_comp[deg_ess_comp$p_val_corr>0.05,"data"],'NA')

ggplot(deg_ess_comp,aes(x= isolate, y =corrected, fill = data_plot))+
  geom_bar(stat = 'identity',position = 'dodge',col='black',alpha=1)+
  theme_classic()+
  facet_grid(~type)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5, 1), labels = c(0,0.5,1,1.5,2))+
  ylab('Odds ratios')+
  xlab('Isolates')+
  scale_fill_manual(values = c('lightsalmon','lightgrey','royalblue','lightgrey'))





####Figure 2####
#####A#####
x= cor(mean_norm_RNA, method = 's')[lower.tri(cor(mean_norm_RNA, method = 's'))]
temp= cor(mean_norm_Ribo, method = 's')[lower.tri(cor(mean_norm_Ribo, method = 's'))]
a = data.frame('RNA-seq'=x, 'Ribo-seq'=temp)
a= melt(a)
ggplot(a,aes(variable, value, fill=variable ))+
  geom_violin(alpha= 0.7)+
  geom_boxplot(fill='white', width= 0.1)+
  theme_classic()+
  stat_compare_means(paired = T,label = 'p.signif',label.x = 1.5)+
  scale_fill_manual(values = c( 'royalblue','lightsalmon'))+
  ylab('Correlation')+
  xlab('Dataset')+
  theme(legend.position = 'none')


wilcox.test(x, temp, paired = T)
median(temp)/median(x)

#####B#####

x= as.matrix(dist(t(mean_norm_RNA)))[lower.tri(as.matrix(dist(t(mean_norm_RNA))))]
temp=  as.matrix(dist(t(mean_norm_Ribo)))[lower.tri(as.matrix(dist(t(mean_norm_Ribo))))]
a = data.frame('RNA-seq'=x, 'Ribo-seq'=temp)
a= melt(a)
ggplot(a,aes(variable, value, fill=variable ))+
  geom_violin()+
  geom_boxplot(fill='white', width= 0.1)+
  theme_classic()+
  stat_compare_means(paired = T,label = 'p.signif',label.x = 1.5)+
  scale_fill_manual(values = c( 'royalblue','lightsalmon'))+
  ylab('Distance')+
  xlab('Dataset')+
  theme(legend.position = 'none')

wilcox.test(x, temp, paired = T)
median(temp)/median(x)

#####C#####

annota_buff_each_strain$CQC@contrasts
anota2seqPlotFC(annota_buff_each_strain$CQC,plotToFile = F,selContrast =4,contrastName = 'CQC vs BPL')
pdf("2C.pdf", width= 5, height = 5) 
# 2. Create a plot
anota2seqPlotFC(annota_buff_each_strain$CQC,plotToFile = F,selContrast =4,contrastName = 'CQC vs BPL')
# Close the pdf file
dev.off()



####Figure 3####
#####A#####

x = stack(comp_new)
colnames(x)=c('ORF', 'Complex')
x= x[x$ORF%in% rownames(norm_Ribo),]

x=x[!x$Complex%in% names(table(x$Complex)[table(x$Complex)==1]),]


temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    c = matrix(ncol=2,nrow=2)
    c[1,1]=length(intersect(unique(unlist(comp_new)),b))
    c[1,2]= length(b)-c[1,1]
    c[2,1]=length(intersect(unique(unlist(comp_new)),d))
    c[2,2]=length(d)-c[2,1]
    rownames(c)=c('transcript_buff','transcript_not_buffered')
    colnames(c)=c('comp', 'not comp')
    c(fisher.test(c)$estimate,fisher.test(c)$p.value)
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

buff_comp= do.call(rbind,temp)

buff_comp$log10_pval=-log10(buff_comp$V2)


buff_comp$sign=buff_comp$V2<0.05
buff_comp$OR_cor= buff_comp$`odds ratio`
buff_comp[!buff_comp$sign,"OR_cor"]=NA
ggplot(buff_comp,aes(strain,vs, col=OR_cor,size= log10_pval))+
  geom_point()+
  theme_test()+
  scale_color_gradient2('Odds ratio',midpoint = 1,low = 'blue',high = 'red',mid = 'white',na.value = 'transparent',lim= c(0,4))+
  scale_size('-Log10(p-value)',breaks = c(0,1,2.5,5,8))+
  ylab('Isolates')+
  xlab('Isolates')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  


#####B#####


temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    c = matrix(ncol=2,nrow=2)
    c[1,1]=length(intersect(unique(ess),b))
    c[1,2]= length(b)-c[1,1]
    c[2,1]=length(intersect(unique(ess),d))
    c[2,2]=length(d)-c[2,1]
    rownames(c)=c('transcript_buff','transcript_not_buffered')
    colnames(c)=c('ess', 'not ess')
    c(fisher.test(c)$estimate,fisher.test(c)$p.value)
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

buff_ess= do.call(rbind,temp)

buff_ess$log10_pval=-log10(buff_ess$V2)


buff_ess$sign=buff_ess$V2<0.05
buff_ess$OR_cor= buff_ess$`odds ratio`
buff_ess[!buff_ess$sign,"OR_cor"]=NA

ggplot(buff_ess,aes(strain,vs, col=OR_cor,size= log10_pval))+
  geom_point()+
  theme_test()+
  scale_color_gradient2('Odds ratio',midpoint = 1,low = 'blue',high = 'red',mid = 'white',na.value = 'transparent',lim= c(0,4))+
  scale_size_continuous('-Log10(p-value)',breaks = c(1,2.5,5,8))+
  ylab('Isolates')+
  xlab('Isolates')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#####C#####

x=rowMeans(norm_RNA)

temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    d=stack(x[d])
    d$type = 'Transcipt variation'
    
    b=stack(x[b])
    b$type = 'Buffered transcipt variation'
    wilcox.test(b$values,d$values)
    c=rbind(b,d)
    c = aggregate(data=c,FUN=median, values~type )
    c(setNames(c$values,c$type),wilcox.test(b$values,d$values)$p.value)
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

buff_expr= do.call(rbind,temp)

buff_expr$log10_pval=-log10(buff_expr$V3)


buff_expr$sign=buff_expr$V3<0.05
temp=(melt(buff_expr[,1:2]))
table(temp$value)

ggplot(temp[!duplicated(temp$value),],aes(variable, value, fill= variable))+
  geom_violin()+
  geom_boxplot(width=0.1, fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4,label = 'p.signif')+
  ylab('Median RNA abundace')+
  xlab('Gene group')+
  scale_fill_manual(values = c('#89ABE3FF','#A8D5BAFF'))+
  theme(legend.position = 'none')


ggplot(temp[!duplicated(temp$value),],aes(variable, value, fill= variable))+
  geom_violin()+
  geom_boxplot(width=0.1, fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4)+
  ylab('Median RNA abundace')+
  xlab('Gene group')+
  scale_fill_manual(values = c('#89ABE3FF','#A8D5BAFF'))+
  theme(legend.position = 'none')
x= temp[!duplicated(temp$value),]
aggregate(data=x, FUN=median, value~variable)
#####D#####


a= mget(ls()[grep('tai2_',ls())])
a=a[-1]

b= lapply(rownames(norm_Ribo),function(i){
  x=c()
  for(j in a){
    x = c(x,j[i])  
  }
  mean(na.omit(x))
})
b=unlist(setNames(b, nm=rownames(norm_Ribo)))

x = b

temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    d=stack(x[d])
    d$type = 'transcipt variation'
    
    b=stack(x[b])
    b$type = 'buffered transcipt variation'
    wilcox.test(b$values,d$values)
    c=rbind(b,d)
    c = aggregate(data=c,FUN=median, values~type )
    c(setNames(c$values,c$type),wilcox.test(b$values,d$values)$p.value)
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

buff_tai= do.call(rbind,temp)

buff_tai$log10_pval=-log10(buff_tai$V3)


buff_tai$sign=buff_tai$V3<0.05
temp=(melt(buff_tai[,1:2]))
table(temp$value)
ggplot(temp[!duplicated(temp$value),],aes(variable, value,fill=variable))+
  geom_violin()+
  geom_boxplot(width=0.1, fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4,label = 'p.signif')+
  ylab('tAI')+
  xlab('Gene group')+
  scale_fill_manual(values = c('#89ABE3FF','#A8D5BAFF'))+
  theme(legend.position = 'none')

aggregate(data=temp[!duplicated(temp$value),],FUN= mean, value~variable)
####Figure 4 ####

######A######
a =introgression_data[,grep('RNA',colnames(introgression_data))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.5, paired = T,label = 'p.signif')+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('ORF')+
  ylab('RNA-seq level')+
  scale_fill_manual(values = c('grey','#00B6EB'))

######B######

a =introgression_data[,grep('Ribo',colnames(introgression_data))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.5, paired = T,label = 'p.signif')+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('ORF')+
  ylab('Ribo-seq level')+
  scale_fill_manual(values = c('grey','#00B6EB'))

#####C#####
orf_info <- fread('Data/orf_info.csv', data.table = F)

orf_info<- orf_info[-166,]
rownames(orf_info)<- as.character(orf_info$`Annotation Name`)
sup_orf = read.csv("data/old_supORF_count/counts.sup_ORF_strains.csv")

j = 1
b = sup_orf[,c(1,j+2)]



rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b
b = rep_data_pan[rep_data_pan$V1%in%x  ,grep('BPL',colnames(rep_data_pan))]
b = b[-which(rowSums(b)==0),]
b=b[!rowSums(b[,1:2])==0,]
b=b[!rowSums(b[,3:4])==0,]
b= tpm_pangenome[rownames(b),]
b[rep_data_pan[rownames(b),colnames(b)]==0]=0


c = rownames(b)[which(orf_info[rownames(b),7]=='HGT')]

b= (b[c,]) 
temp= b[,grep('BPL',colnames(b))]
temp = rbind(temp,cbind(norm_RNA,norm_Ribo)[,grep('BPL',colnames(cbind(norm_RNA,norm_Ribo)))])
temp=as.data.frame(temp)
temp$HGT = rownames(temp)%in%rownames(b)
temp$TE_1=temp$Ribo_BPL_1/temp$RNA_BPL_1

temp$TE_2=temp$Ribo_BPL_2/temp$RNA_BPL_2

temp$TE_mean=rowMeans(temp[,c("TE_1","TE_2")])
temp$Type = ifelse(temp$HGT,'HGT','Canonical ORF')
ggplot(temp,aes(Type, TE_mean, fill=Type))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4,label = 'p.signif')+
  scale_y_log10()+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c('grey','#FAA094FF'))+
  xlab('ORF')+
  ylab('Average translation efficiency')

table(temp$Type)
aggregate(data=temp, FUN = mean, TE_mean~Type)
####_####
####Sup figures ####
####_####

####Figure S2####
#####A#####
x= cor(norm_RNA, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  theme_test()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value), lim = c(0.4,1), name='Correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,3)), size = 4)+
  xlab('')+
  ylab('')

#####B#####
x= cor(norm_Ribo, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  theme_test()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value), lim = c(0.4,1),name='Correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,3)), size = 4)+
  xlab('')+
  ylab('')

####Figure S3####
#####A#####
x=t(norm_RNA)
x=as.data.frame(x)
x$strain = substr(rownames(x),5,7 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,label = F)+
  theme_classic()

#####B#####

x=t(norm_Ribo)
x=as.data.frame(x)
x$strain = substr(rownames(x),6,8 )
autoplot(prcomp(x[, -((ncol(x)-1):ncol(x))]), data = x, col = 'strain' ,label = F)+
  theme_classic()

####Figure S4####
a = lapply(isolates, function(i){
  c(cor.test(mean_norm_RNA[,i], mean_norm_Ribo[,i],method = 's')$estimate,i)
})

a=do.call(rbind,a)
a= as.data.frame(a)
a$rho=as.numeric(a$rho)
ggplot(a,aes(V2,rho, fill=V2))+
  geom_bar(stat = 'identity',alpha = 0.7)+
  xlab('Isolate')+
  ylab('Rho')+
  theme_classic()+
  theme(legend.position = 'none')

####Figure S5####
#####A#####
fit <- lmFit(norm_Ribo, design)
CUTOFF_LFC = log2(1.2)
one_contrast <- makeContrasts(contrasts = df_group$strain,levels = design)
one_contrast[one_contrast<1] <- -1 / (nlevels(df_group$strain)-1) # weighted mean expression from all strains
one_fit <- contrasts.fit(fit, one_contrast)
one_DE <- treat(one_fit,trend=T,robust=T,fc=1.2) # In practice, significantly DE gene will have a fc much higher than 2

CUTOFF_LFC = log2(1.2)

temp <- pblapply(df_strains$strain, function(i){ref_strain=i


strain_one_DE = topTreat(fit = one_DE,  coef = ref_strain , lfc=0, p.value = 1, adjust.method = 'fdr', number = Inf) |> 
  mutate(DE = case_when(logFC < (-1*CUTOFF_LFC) ~ "down",  logFC > CUTOFF_LFC ~ "up"))
strain_one_DE$isolate=i
strain_one_DE$gene = rownames(strain_one_DE)
return(strain_one_DE)
})
names(temp)= df_strains$strain
one_vs_all_Ribo=do.call(rbind,temp)
one_vs_all_Ribo$is_DE = one_vs_all_Ribo$adj.P.Val<1e-5
sum(one_vs_all_Ribo$is_DE)
table(one_vs_all_Ribo[one_vs_all_Ribo$is_DE,"isolate"])


DE_gene = unique(one_vs_all_Ribo[one_vs_all_Ribo$is_DE,"gene"])

c = list()

x= pblapply(df_strains$strain, function(i){
  
  a=temp[[i]]
  a$DEP = a$adj.P.Val < 1e-5
  p = ggplot(a,aes(logFC, -log10(P.Value),col = DEP))+
    geom_point()+
    scale_color_manual(values = c('black','purple1'))+
    theme_classic()+
    theme(legend.position = 'none')+
    ggtitle(i)
  #plot(p)
  c<<-list.append(c,p)
  a$strain = i
  return(a)
})
do.call(grid.arrange, c(c,ncol=4))


#####B#####


fit <- lmFit(norm_RNA, design)
CUTOFF_LFC = log2(1.2)
one_contrast <- makeContrasts(contrasts = df_group$strain,levels = design)
one_contrast[one_contrast<1] <- -1 / (nlevels(df_group$strain)-1) # weighted mean expression from all strains
one_fit <- contrasts.fit(fit, one_contrast)
one_DE <- treat(one_fit,trend=T,robust=T,fc=1.2) # In practice, significantly DE gene will have a fc much higher than 2

CUTOFF_LFC = log2(1.2)

temp <- pblapply(df_strains$strain, function(i){ref_strain=i


strain_one_DE = topTreat(fit = one_DE,  coef = ref_strain , lfc=0, p.value = 1, adjust.method = 'fdr', number = Inf) |> 
  mutate(DE = case_when(logFC < (-1*CUTOFF_LFC) ~ "down",  logFC > CUTOFF_LFC ~ "up"))
strain_one_DE$isolate=i
strain_one_DE$gene = rownames(strain_one_DE)
return(strain_one_DE)
})
names(temp)= df_strains$strain
one_vs_all_RNA=do.call(rbind,temp)
one_vs_all_RNA$is_DE = one_vs_all_RNA$adj.P.Val<1e-5
sum(one_vs_all_RNA$is_DE)
table(one_vs_all_RNA[one_vs_all_RNA$is_DE,"isolate"])


DE_gene = unique(one_vs_all_RNA[one_vs_all_RNA$is_DE,"gene"])

c = list()

x= pblapply(df_strains$strain, function(i){
  
  a=temp[[i]]
  a$DEP = a$adj.P.Val < 1e-5
  p = ggplot(a,aes(logFC, -log10(P.Value),col = DEP))+
    geom_point()+
    scale_color_manual(values = c('black','purple1'))+
    theme_classic()+
    theme(legend.position = 'none')+
    ggtitle(i)
  #plot(p)
  c<<-list.append(c,p)
  a$strain = i
  return(a)
})
do.call(grid.arrange, c(c,ncol=4))


####Figure S6####
ggplot(deg_ess_comp, aes(y= OR, x= data))+
  geom_boxplot(fill='grey')+
  stat_compare_means(label='p.signif', label.x = 1.5,paired = T)+
  theme_classic()
ggplot(deg_ess_comp, aes(y= OR, x= data))+
  geom_boxplot(fill='grey')+
  stat_compare_means( label.x = 1.5,paired = T)+
  theme_classic()

####Figure S7####

a = data.frame('RNA-seq'=melt(abs_log_RNA)$value, 'Ribo-seq'=melt(abs_log_Ribo)$value)
a= melt(a)
ggplot(a,aes(variable, value, fill=variable ))+
  geom_violin(alpha= 0.7)+
  geom_boxplot(fill='white', width= 0.1)+
  theme_classic()+
  stat_compare_means(label = 'p.signif',label.x = 1.5)+
  scale_fill_manual(values = c( 'royalblue','lightsalmon'))+
  ylab('|log2(FC)|')+
  xlab('Dataset')+
  theme(legend.position = 'none')+
  scale_y_log10()


wilcox.test(melt(abs_log_RNA)$value, melt(abs_log_Ribo)$value)
mean(na.omit(melt(abs_log_RNA)$value))
mean(melt(abs_log_Ribo)$value)

#### Figure S8####
#####A#####
x = pblapply(isolates, function(i){
  x = annota_buff_each_strain[[i]]@selectedBuffering@selectedRvmData
  temp = lapply(x, function(temp){
    temp<<-temp
    rownames(temp[temp$singleRegMode=='buffering',])
  })
})

names(x)=isolates

buffered_genes_anota = x


temp= stack(lapply(buffered_genes_anota, function(i){i<<-i
lapply(i, length)
}))
temp$ind=as.character(temp$ind)
ggplot(temp,aes(ind, values, fill=ind))+
  geom_boxplot()+
  ylab('Number of buffered genes')+
  xlab('Isolates')+
  theme_classic()+
  theme(legend.position = 'none')
range(temp$values)
mean(temp$values)


#####B#####
temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    length(b)/(length(b)+length(d))
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

temp=do.call(rbind,temp)
ggplot(temp,aes(strain, V1,fill=strain))+
  geom_boxplot()+
  ylab('Proportion of buffered genes')+
  xlab('Isolates')+
  theme_classic()+
  theme(legend.position = 'none')
####Figure S9####

x = buff_comp[!duplicated(buff_comp$`odds ratio`),c(1,2)]
x$type ='Complex'
temp = x

x = buff_ess[!duplicated(buff_ess$`odds ratio`),c(1,2)]
x$type ='Essential'
temp = rbind(x,temp)
ggplot(temp, aes(y=type, x = `odds ratio`))+
  geom_boxplot(fill='lightgrey',outlier.colour = 'transparent')+
  geom_jitter(height = 0.1)+
  theme_classic()+
  scale_x_continuous(breaks = c(0,1,2,3,4), limits = c(0,4))+
  geom_vline(xintercept = 1)+
  ylab('Gene group')+
  xlab('Odds ratios')
  

wilcox.test(temp[temp$type=='Complex',"odds ratio"],mu = 1)
wilcox.test(temp[temp$type=='Essential',"odds ratio"],mu = 1)


####Figure S10####
x=rowMeans(norm_RNA)

temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    d=stack(x[d])
    d$type = 'transcipt variation'
    
    b=stack(x[b])
    b$type = 'buffered transcipt variation'
    wilcox.test(b$values,d$values)
    c=rbind(b,d)
    
    c$p_value=wilcox.test(b$values,d$values)$p.value
    c$vs = colnames(a@contrasts)[i]
    c
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b
})

temp= do.call(rbind,temp)
temp$label = paste(temp$type,temp$vs,temp$strain)

temp$comparison = paste(temp$strain,temp$vs,sep = ' vs ')
x = combn(isolates,2)
x=t(x)
x=as.data.frame(x)
x = paste(x$V1,x$V2,sep = ' vs ')
temp=temp[temp$comparison%in%x,]
ggplot(temp, aes(x=type, fill = type, y = values))+
  geom_violin()+
  geom_boxplot(fill='white', width=0.1)+
  theme_classic()+
  facet_wrap(~comparison,ncol = 7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(label = 'p.signif',label.y = 10, label.x = 1.5)+
  scale_fill_manual('Gene group',values = c('#89ABE3FF','#A8D5BAFF'))+
  stat_summary(fun=median, geom="line", aes(group=1))+
  xlab('Gene group')+
  theme(axis.text.x = element_blank())+
  ylab('mRNA level')


####Figure S11####

x=rowMeans(norm_RNA)

temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = rownames(norm_Ribo)
    d=d[!d%in%b]
    d=stack(x[d])
    d$type = 'Other genes'
    
    b=stack(x[b])
    b$type = 'Buffered transcipt variation'
    wilcox.test(b$values,d$values)
    c=rbind(b,d)
    c = aggregate(data=c,FUN=median, values~type )
    c(setNames(c$values,c$type),wilcox.test(b$values,d$values)$p.value)
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

temp= do.call(rbind,temp)

temp$log10_pval=-log10(temp$V3)


temp$sign=temp$V3<0.05
temp=(melt(temp[,1:2]))
table(temp$value)

ggplot(temp[!duplicated(temp$value),],aes(variable, value, fill= variable))+
  geom_violin()+
  geom_boxplot(width=0.1, fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4,label = 'p.signif')+
  ylab('Median RNA abundace')+
  xlab('Gene group')+
  scale_fill_manual(values = c('#89ABE3FF','grey'))+
  theme(legend.position = 'none')

####Figure S12####

x = data.frame(mRNA.levels=rowMeans(norm_RNA), mRNA.variance = apply(norm_RNA,1, var))
ggplot(x, aes(`mRNA.levels`,`mRNA.variance`))+
  geom_point()+
  theme_classic()+
  stat_cor(method = 's')+
  scale_y_log10()+
  xlab('mRNA levels')+
  ylab('mRNA variance')

####Figure S13####

a= mget(ls()[grep('tai2_',ls())])

b= lapply(rownames(norm_Ribo),function(i){
  x=c()
  for(j in a){
    x = c(x,j[i])  
  }
  mean(na.omit(x))
})
b=unlist(setNames(b, nm=rownames(norm_Ribo)))

x = b


temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='buffering',]
    b = rownames(b)
    d = a@selectedTotalmRNA@selectedRvmData[[i]]
    d= d[d$singleRegMode=='abundance',]
    d=rownames(d)
    d=stack(x[d])
    d$type = 'transcipt variation'
    
    b=stack(x[b])
    b$type = 'buffered transcipt variation'
    wilcox.test(b$values,d$values)
    c=rbind(b,d)
    
    c$p_value=wilcox.test(b$values,d$values)$p.value
    c$vs = colnames(a@contrasts)[i]
    c
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b
})

temp= do.call(rbind,temp)
temp$label = paste(temp$type,temp$vs,temp$strain)

temp$comparison = paste(temp$strain,temp$vs,sep = ' vs ')
x = combn(isolates,2)
x=t(x)
x=as.data.frame(x)
x = paste(x$V1,x$V2,sep = ' vs ')
temp=temp[temp$comparison%in%x,]
ggplot(temp, aes(x=type, fill = type, y = values))+
  geom_violin()+
  geom_boxplot(fill='white', width=0.1)+
  theme_classic()+
  facet_wrap(~comparison,ncol = 7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(label = 'p.signif', label.x = 1.5,label.y = 0.5)+
  scale_fill_manual('Gene group',values = c('#89ABE3FF','#A8D5BAFF'))+
  stat_summary(fun=median, geom="line", aes(group=1))+
  xlab('Gene group')+
  theme(axis.text.x = element_blank())+
  ylab('tAI')

####Figure S14####



d = unique(orf_info$`Origin assignment`)
d = data.frame( gg_color_hue(7),row.names = d)
ori <-d

e = matrix(ncol = 3)
e=e[-1,]

for(i in 1:8){
  
  b = sup_orf[,c(1,i+2)]
  rownames(b)<- b$X
  b <-rownames(b)[which(!is.na(b[,2]))]
  
  b = tpm_pangenome[b,grep(isolates[i],colnames(tpm_pangenome))]
  
  a = (table(orf_info[rownames(b),7]))
  a = as.data.frame(a)
  a$Expr = 'All'
  
  d<-a
  
  
  
  b = sup_orf[,c(1,i+2)]
  rownames(b)<- b$X
  b <-rownames(b)[which(!is.na(b[,2]))]
  b = tpm_pangenome[b,grep(isolates[i],colnames(tpm_pangenome))]
  b[rep_data_pan[rownames(b),colnames(b)]==0]=0
  b = b[!(rowSums(b)==0),]
  b=b[!b[,1]==0,]
  b=b[!b[,2]==0,]
  a = (table(orf_info[rownames(b),7]))
  a = as.data.frame(a)
  a$Expr = 'Expressed'
  
  b = rbind(d,a)
  
  b = b[order(b$Var1),]
  b= b[order(b$Expr),]
  d=spread(b, Expr, Freq)
  d[is.na(d)]=0
  d$All=cumsum(rev(d$All))
  d$Expressed=cumsum(rev(d$Expressed))
  c = fct_inorder(b$Var1)
  c =as.character(c)
  c=as.character( ori[c,])
  
  e = rbind(e,c(isolates[i],sum(b[b$Expr=='All',"Freq"]),sum(b[b$Expr=='Expressed',"Freq"])))
  
  p = ggplot() +
    geom_bar(data = b,
             aes(x = Expr, y =Freq, fill = Var1),
             colour = 'black', width = 0.3, stat="identity") +
    geom_segment(data = d,
                 colour = "black",
                 aes(x = 1 + 0.3/2,
                     xend = 2 - 0.3/2,
                     y = All,
                     yend = Expressed)) +
    theme(panel.background = element_rect(fill = "white"),        
          panel.grid = element_blank()) +
    scale_fill_manual(values=c)+
    theme_classic()+
    theme(legend.position="none")+
    ggtitle(isolates[i])+
    labs(fill = "  ", y=NULL, x= NULL)
  assign(paste('p',i,sep = ''),p)
  print(e) 
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4 )


x= ori
x$ori=rownames(x)

temp <- cowplot::get_legend(ggplot(x, aes(x=1, y=ori, fill=factor(ori, levels = ori)))+
                              geom_tile()+
                              scale_fill_manual('ORF origin',values=ori$gg_color_hue.7.))

grid.newpage()
grid.draw(temp)

####Figure S15####
#####A#####
orf_info <- fread('Data/orf_info.csv', data.table = F)

orf_info<- orf_info[-166,]
rownames(orf_info)<- as.character(orf_info$`Annotation Name`)
sup_orf = read.csv("data/old_supORF_count/counts.sup_ORF_strains.csv")

j = 1
b = sup_orf[,c(1,j+2)]

rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b
b = rep_data_pan[rep_data_pan$V1%in%x  ,grep('BPL',colnames(rep_data_pan))]
b = b[-which(rowSums(b)==0),]
b=b[!rowSums(b[,1:2])==0,]
b=b[!rowSums(b[,3:4])==0,]
b= tpm_pangenome[rownames(b),]
b[rep_data_pan[rownames(b),colnames(b)]==0]=0


c = rownames(b)[which(orf_info[rownames(b),7]=='HGT')]
b= (b[c,]) 
temp= b[,grep('BPL',colnames(b))]
temp = rbind(temp,cbind(norm_RNA,norm_Ribo)[,grep('BPL',colnames(cbind(norm_RNA,norm_Ribo)))])
temp=as.data.frame(temp)
temp$HGT = rownames(temp)%in%rownames(b)
temp$ribo = rowMeans(temp[,c(3,4)])
temp$RNA = rowMeans(temp[,c(1,2)])

ggplot(temp[!temp$HGT,],aes(RNA,ribo))+
  geom_point()+
  theme_classic()+
  geom_point(data=temp[temp$HGT,],aes(RNA,ribo), col='salmon')+
  stat_cor(method = 's')+
  geom_abline(intercept =  0,slope =1)+
  xlab('mRNA level')+
  ylab('Ribo-seq level')

#####B#####

orf_info <- fread('Data/orf_info.csv', data.table = F)

orf_info<- orf_info[-166,]
rownames(orf_info)<- as.character(orf_info$`Annotation Name`)
sup_orf = read.csv("data/old_supORF_count/counts.sup_ORF_strains.csv")

j = 1
b = sup_orf[,c(1,j+2)]



rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b

for(i in 'BPL'){
  
  require("tAI")
  eco.trna <- get(paste(i, '.trna2', sep = ''))
  eco.ws <- get.ws(tRNA=eco.trna, sking=1)
  eco.m <- matrix(scan('Data/8_cds/misc_rep/allORFs_pangenome.m'), ncol=61, byrow=TRUE)
  eco.m <- eco.m[,-33]
  eco.tai <- get.tai(eco.m, eco.ws)
  temp =read.fasta("Data/8_cds/misc_rep/allORFs_pangenome.ffn")
  
  temp = getName(temp)
  temp = setNames( eco.tai, temp)
  assign(paste('tai2_allORF_',i,sep=''), temp)
}

orf_info <- fread('Data/orf_info.csv', data.table = F)

orf_info<- orf_info[-166,]
rownames(orf_info)<- as.character(orf_info$`Annotation Name`)
sup_orf = read.csv("data/old_supORF_count/counts.sup_ORF_strains.csv")

j = 1
b = sup_orf[,c(1,j+2)]



rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b
b = rep_data_pan[rep_data_pan$V1%in%x  ,grep('BPL',colnames(rep_data_pan))]
b = b[-which(rowSums(b)==0),]
b=b[!rowSums(b[,1:2])==0,]
b=b[!rowSums(b[,3:4])==0,]
b= tpm_pangenome[rownames(b),]
b[rep_data_pan[rownames(b),colnames(b)]==0]=0


c = rownames(b)[which(orf_info[rownames(b),7]=='HGT')]

b= (b[c,]) 
temp= b[,grep('BPL',colnames(b))]
temp = rbind(temp,cbind(norm_RNA,norm_Ribo)[,grep('BPL',colnames(cbind(norm_RNA,norm_Ribo)))])
temp=as.data.frame(temp)
temp$HGT = rownames(temp)%in%rownames(b)
temp$TE_1=temp$Ribo_BPL_1/temp$RNA_BPL_1

temp$TE_2=temp$Ribo_BPL_2/temp$RNA_BPL_2

temp$TE_mean=rowMeans(temp[,c("TE_1","TE_2")])
temp$Type = ifelse(temp$HGT,'HGT','Canonical ORF')

x= rbind(data.frame(tAI=tai2_allORF_BPL[rownames(temp)[temp$Type=='HGT']], type='HGT'),data.frame(tAI=tai2_BPL, type='cannonical'))
ggplot(x,aes(type,tAI,fill=type))+
  geom_violin()+
  geom_boxplot(width=0.1,fill='white')+
  theme_classic()+
  stat_compare_means(label = 'p.signif', label.x = 1.5)+
  scale_fill_manual('Gene group',values = c('grey','salmon'))+
  theme(legend.position = 'none')+
  xlab('ORF')

####_####
####Sup table ####
####_####
####Table S3####
write_clip(norm_RNA)

####Table S4####
write_clip(norm_Ribo)
####Table S5####
write_clip(norm_RNA_pan)

####Table S6####
write_clip(norm_Ribo_pan)

####Table S7####
x=one_vs_all_RNA
x$data='RNA-seq'
temp=one_vs_all_Ribo
temp$data='Ribo-seq'
x = rbind(x,temp)
write_clip(x)


####Table S8####
gostres <- gost(query = unique(one_vs_all_Ribo[one_vs_all_Ribo$is_DE,"gene"]) , organism = "scerevisiae", 
                ordered_query = FALSE,  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, user_threshold = 0.05, correction_method = "g_SCS",  
                domain_scope = "custom", custom_bg =rownames( norm_RNA), numeric_ns = "", sources = NULL, as_short_link = FALSE) 
a = (gostres$result)
a = a[grep('GO',a$source),]
temp = a
gostres <- gost(query = unique(one_vs_all_RNA[one_vs_all_RNA$is_DE,"gene"]) , organism = "scerevisiae", 
                ordered_query = FALSE,  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, user_threshold = 0.05, correction_method = "g_SCS",  
                domain_scope = "custom", custom_bg =rownames( norm_RNA), numeric_ns = "", sources = NULL, as_short_link = FALSE) 
a = (gostres$result)
a = a[grep('GO',a$source),]

temp$query='Ribo-seq'
a$query='RNA-seq'
x = rbind(a,temp)
write_clip(x[,-14])

#### Table S9 ####


temp= lapply(isolates,function(a){
  a<<-a
  x= a
  a= annota_buff_each_strain[[a]]
  b= lapply(1:7,function(i){
    b = a@selectedBuffering@selectedRvmData[[i]]
    b$contrast = paste0(x,'vs',colnames(a@contrasts)[i])
    b=b[b$singleRegMode!='translation',]
    b
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  
})
x = do.call(rbind,temp)
write_clip(x)

####Table S10####

temp= lapply(annota_buff_each_strain,function(a){
  a<<-a
  
  b= lapply(1:7,function(i){
    b = a@selectedTranslation@selectedRvmData[[i]]
    
    b = b[b$singleRegMode=='translation',]
    b = rownames(b)
    
    c = a@selectedBuffering@selectedRvmData[[i]]
    
    c = c[c$singleRegMode=='buffering',]
    c = rownames(c)
    
    e = a@selectedTotalmRNA@selectedRvmData[[i]]
    
    e = e[e$singleRegMode=='abundance',]
    e = rownames(e)
    
    
    c(length(c),length(b),length(e))
  })
  b=do.call(rbind,b)
  b=as.data.frame(b)
  b$strain= isolates[!isolates%in%colnames(a@contrasts)]
  b$vs = colnames(a@contrasts)
  b
})

a= do.call(rbind, temp)
a$unique = !duplicated(paste(a$V1,a$V2))

a=a[a$unique,]
a$ratio = a$V1/a$V2
mean(a$ratio)
a$prop = a$V1/(a$V1+a$V3)
colnames(a)= c('Number of buffered genes',     "Number of enhanced genes"  , "Number of genes with according transcription and translation variation" ,  "strain" ,"vs"   ,  "unique comparison", "buffering/enhanced ratio",'Proportion of buffered genes among genes with transcriptional variation' )
write_clip(a)
####Table S11####


temp= stack(lapply(buffered_genes_anota, function(i){i<<-i
unlist(i)
}))

x = as.data.frame(table(temp$values)/2)
x$Var1= as.character(x$Var1)
x=rbind(x,data.frame(Var1=rownames(norm_Ribo)[!rownames(norm_Ribo)%in%x$Var1],Freq= rep(0,length(rownames(norm_Ribo)[!rownames(norm_Ribo)%in%x$Var1]))))

temp=setNames(x$Freq,nm = x$Var1)
temp=sort(temp)
temp= fgsea(temp,pathways=GO2geneID )

temp = na.omit(temp)
temp= temp[temp$padj<0.1,]
temp=temp[order(temp$pathway),]
temp$term=go2term(temp$pathway)[,2]
write_clip(temp[order(temp$NES),-c(8)])


####Table S12####


a= mget(ls()[grep('tai2_',ls())])
x= lapply(isolates, function(i){
  b=grep(i,names(a))
  b=a[[b]]
  c=intersect(names(b),rownames(norm_Ribo))
  c = data.frame(mRNA_level= mean_norm_RNA[c,i],tAI=b[c])
  c(as.numeric(cor.test(c$mRNA_level,c$tAI,method = 's')$estimate),i,'RNA-seq', cor.test(c$mRNA_level,c$tAI,method = 's')$p.value)
})

temp=as.data.frame(do.call(rbind,x))

a= mget(ls()[grep('tai2_',ls())])
x=lapply(isolates, function(i){
  b=grep(i,names(a))
  b=a[[b]]
  c=intersect(names(b),rownames(norm_Ribo))
  c = data.frame(ribo_level= mean_norm_Ribo[c,i],tAI=b[c])
  c(as.numeric(cor.test(c$ribo_level,c$tAI,method = 's')$estimate),i,'Ribo-seq', cor.test(c$ribo_level,c$tAI,method = 's')$p.value)
})
x=rbind(temp,as.data.frame(do.call(rbind,x)))
colnames(x)=c('Correlation coefficient', 'isolate','data','pvalue')
write_clip(x)
