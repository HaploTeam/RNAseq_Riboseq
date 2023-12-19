####Package####

{ 
  library(DESeq2)
  library(stringr)
  library(gprofiler2)
  library(FSA)
  library(limma)
  library(corrplot)
  library(ggVennDiagram)
  library(Rtsne)
  library(fgsea)
  library(KEGGREST)
  library(org.Sc.sgd.db)
  library(nortest)
  library(rrvgo)
  #library('ROTS')
  library(impute)
  library(mgcv)
  library(ggtree)
  library(pegas)
  #library(edgeR)
  library(magick)
  library(affy)
  library(relaimpo)
  library(miscTools)
  library(fitdistrplus)
  library('MALDIquant')
  library(multcompView)
  library(topGO)
  library(gtools)
  library(pbapply)
  #library(WGCNA)
  library(Biostrings)
  library(seqinr)
  library(intrval)
  library(gridExtra)
  library(rlist)
  library(zoom)
  library(GGally)
  library(dendextend)
  library(factoextra)
  library(NbClust)
  library(clipr)
  library(dplyr)
  library(bigutilsr)
  library(performance) 
  library(ggplot2)
  library(ggfortify)
  library(tidyverse)
  library(itertools)
  library(gplots)
  library(ape)
  library(clusterProfiler)
  library(enrichplot)
  library(data.table)
  library(ggpubr)
  library(BiocManager)
  library(readr)
  library(ggcorrplot)
  library(plotly)
  library(ggsignif)
  #library(pathfindR)
  library(anota2seq)
 
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  

  
}
####data ####
gene_info = fread('Data/CompleteGeneAnnot_strand.csv', data.table = F)
rownames(gene_info)<-gene_info$systematic_name
download_Yeast_GO_mapping <-
  function(yeast.GO.url='http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'){
    ##### download Yeast Gene Ontology mapping  http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz
    GO.assocs.all  <- read.table(textConnection(readLines(gzcon(url(yeast.GO.url))))
                                 ,header=FALSE,check.names=FALSE, stringsAsFactors=FALSE, sep="\t", quote=NULL, comment.char ='!'
    )
    my.col.names <- c("Database","SGDID","DB_Object_Symbol","NOT","GOID","DB:Ref","Evidence"
                      ,"With:From","GO_Aspect","DB_Object_Name","Synonym","Type","Taxon","Date","Asigned By","Notes 1")
    GO.assocs.all <- GO.assocs.all[,1:length(my.col.names)] ### sometimes there are empty extra columns
    colnames(GO.assocs.all) <- my.col.names
    
    ### use only most important columns "SGDID", "DB_Object_Symbol", "GOID"
    GO.assocs <- GO.assocs.all[,match(c("SGDID", "GOID"),colnames(GO.assocs.all))]
    ### one can filter for a certain evidence like IDA : Inferred from Direct Assay 
    ### in our case we will accept any of the Evidences and keep unique records of the association
    GO.assocs  <- unique(GO.assocs) #### for some reasons we have same rows duplications (for example the multiple Evidence levels: IBA, IC, IDA etc)
  }

Yeast.GO.assocs =download_Yeast_GO_mapping(  yeast.GO.url = "http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz")

Yeast.GO.assocs=Yeast.GO.assocs[Yeast.GO.assocs$SGDID%in%gene_info$SGD_name,]

a = unique(Yeast.GO.assocs$SGDID)
b = lapply(a,function(i){
  b = which(gene_info$SGD_name%in%i)
  return(b)
})
b= as.character(b)
b = rownames(gene_info)[as.numeric(b)]
a= cbind(a,b)
a= as.data.frame(a)
rownames(a)=a$a
Yeast.GO.assocs$STDID= a[Yeast.GO.assocs$SGDID,"b"]
a =  unique(Yeast.GO.assocs$STDID)
b = lapply(a, function(i){
  b = Yeast.GO.assocs[Yeast.GO.assocs$STDID==i,"GOID"]
})
names(b)= unique(Yeast.GO.assocs$STDID)
geneID2GO <- b   

prot_complex <- readxl::read_xls('data/CYC2008_complex.xls')
prot_complex <- as.data.frame(prot_complex)
comp <- unique(prot_complex[,1])
ess_genes_FY <- read.csv('data/reshuffleAnnot.csv', sep = ";")
rownames(ess_genes_FY) <- ess_genes_FY[,1]
ess_genes_FY= na.omit(ess_genes_FY)
ess<- ess_genes_FY[(ess_genes_FY$AnnotFinal=='Yes'),1]
GO2geneID = pblapply(unique(Yeast.GO.assocs$GOID),function(i){
  return(Yeast.GO.assocs[Yeast.GO.assocs$GOID==i,"STDID"])
})
names(GO2geneID)= unique(Yeast.GO.assocs$GOID)


new_complex= fread('Data/complex_costanzo_2016.csv',data.table = F)
comp_new = new_complex$`ORFs annotated to complex`
comp_new= lapply(comp_new,function(i){
  i= str_replace_all(i,' ','')
  i = strsplit(i,';')
  i= unlist(i)
})
names(comp_new)= new_complex$`Protein Complex Name`



####Compute tAI####


{
  for(j in list.files('Data/assembled/tab/csv/')){
    con <- read.csv(paste('Data/assembled/tab/csv/',j, sep=''), sep = '\t')
    colnames(con)[6]= 'anti codon'
    con$`anti codon`=as.character(con$`anti codon`)
    con$type = NA
    
    for(i in 1:nrow(con)){
      a = con[i,6]
      a = DNAString(a)
      a = reverseComplement(a)
      a = as.character(a)
      gsub('T','U',a)
      con[i, 'type'] = a
      
    }
    assign(paste('new_tRNA_CNV', substr(j,1,3), sep='_'),con)
  }
  
  {j= 'S288c.csv'
    con <- read.csv(paste('Data/assembled/S288C/ncbi-genomes-2020-06-12/',j, sep=''), sep = '\t')
    colnames(con)[6]= 'anti codon'
    con$`anti codon`=as.character(con$`anti codon`)
    con$type = NA
    
    for(i in 1:nrow(con)){
      a = con[i,6]
      a = DNAString(a)
      a = reverseComplement(a)
      a = as.character(a)
      gsub('T','U',a)
      con[i, 'type'] = a
      
    }
    assign(paste('new_tRNA_CNV', substr(j,1,3), sep='_'),con)
  }
  
}







a = expand.grid(rep(list(c('T','C','A','G')), 3))
a[,1] = as.character(a[,1])
a[,2]= as.character(a[,2])
a[,3]= as.character(a[,3])
b = c()
for(i in 1:nrow(a)){
  c = paste(a[i,3],a[i,2],a[i,1], sep =  '')
  b=c(b,c)
}

codon_list = b
for(i in  isolates){
  g = get(paste('new_tRNA_CNV_',i, sep=''))
  codon  = unique(g$type)
  codon = na.omit(codon)
  g= na.omit(g)
  temp2=c() 
  for(k in  codon){
    a = sum(g$type==k)
    temp2 = c(temp2, setNames(a, k))
  }
  temp3= toupper( names(temp2))
  for(l in 1:3){temp3=sub("U",'T',temp3)}
  names(temp2) <- temp3
  temp = c()
  for(j in codon_list){
    if(j %in% names(temp2)){temp= c(temp, as.numeric(temp2[j]))}
    if(j %in% names(temp2) ==F){temp= c(temp, 0)}
  }
  assign(paste(i,'.trna2',sep = ''), temp)
}

for(i in isolates){
  
  require("tAI")
  eco.trna <- get(paste(i, '.trna2', sep = ''))
  eco.ws <- get.ws(tRNA=eco.trna, sking=1)
  eco.m <- matrix(scan(paste('Data/8_cds/misc_rep/',i,'.m',sep = '')), ncol=61, byrow=TRUE)
  eco.m <- eco.m[,-33]
  eco.tai <- get.tai(eco.m, eco.ws)
  temp =read.fasta(paste("Data/8_cds/misc_rep/",i,'.ffn',sep = ''))
  
  temp = getName(temp)
  temp = setNames( eco.tai, do.call(rbind,strsplit(temp,'\\|'))[,3])
  assign(paste('tai2_',i,sep=''), temp)
}

remove(temp2)
remove(temp3)


#### QC or bias####
x= cor(ribo_seq, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)

x= cor(RNA_seq, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)

x= cor(norm_RNA, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)

x= cor(norm_Ribo, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)

x=t(norm_Ribo)
x=as.data.frame(x)
x$strain = substr(rownames(x),6,8 )
x$test = rownames(x)
autoplot(prcomp(x[, -((ncol(x)-1):ncol(x))]), data = x, col = 'strain' ,label = TRUE)+
  theme_classic()
  
x=t(norm_RNA)
x=as.data.frame(x)
x$strain = substr(rownames(x),5,7 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,label = TRUE)+
  theme_classic()

x=t(RNA_seq)
x=as.data.frame(x)
x$strain = substr(rownames(x),5,7 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,label = TRUE)

x=t(ribo_seq)
x=as.data.frame(x)
x$strain = substr(rownames(x),6,8 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,label = TRUE)



a = length(intersect(ess,rownames(norm_Ribo)))
b = nrow(norm_Ribo)-a
c =  length(intersect(ess,rownames(ribo_seq)))-a
d = nrow(ribo_seq)-a-b-c
fisher.test(matrix(c(a,b,c,d), nrow = 2))
a = length(intersect(comp,rownames(norm_Ribo)))
b = nrow(norm_Ribo)-a
c =  length(intersect(comp,rownames(ribo_seq)))-a
d = nrow(ribo_seq)-a-b-c
fisher.test(matrix(c(a,b,c,d), nrow = 2))

gostres <- gost(query = rownames(norm_RNA) , organism = "scerevisiae", 
                ordered_query = FALSE,  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, user_threshold = 0.05, correction_method = "g_SCS",  
                domain_scope = "custom", custom_bg =rownames(ribo_seq), numeric_ns = "", sources = NULL, as_short_link = FALSE) 
gostres=gostres$result
gostres=gostres[gostres$source=='GO:BP',]
gostres=gostres[-grep('metabolic', gostres$term_name),]




x= cor(mean_norm_RNA, method = 's')[lower.tri(cor(mean_norm_RNA, method = 's'))]
temp= cor(mean_norm_Ribo, method = 's')[lower.tri(cor(mean_norm_Ribo, method = 's'))]
plot(temp, x, xlab='Ribo_pairwise_cor',ylab='RNA_pairwise_cor')


#####across gene correlation#####

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
range(a$rho)  
#### old data comparison ####


rep_data = fread("canonical_mapping/csv/counts.csv", data.table = F, header = T)
RNA_seq = rep_data[,grep('RNA',colnames(rep_data))]
ribo_seq = rep_data[,grep('Ribo',colnames(rep_data))]
rownames(RNA_seq)= rep_data$V1
rownames(ribo_seq)= rep_data$V1
unique(substr(colnames(ribo_seq),6,8))==unique(substr(colnames(RNA_seq),5,7))
isolates = unique(substr(colnames(RNA_seq),5,7))


old_data_TPM = cbind(fread('/Users/elieteyssonniere/Seafile/dat/dat/Riboseq_8strains/Ribo.TPM.genes.csv', data.table = F),fread('/Users/elieteyssonniere/Seafile/dat/dat/Riboseq_8strains/RNA.TPM.genes.csv', data.table = F))
rownames(old_data_TPM)=old_data_TPM$V1
old_data_TPM=old_data_TPM[-c(1,10:12,21,22)]
old_Ribo_TPM = old_data_TPM[, grep('Ribo', colnames(old_data_TPM))]
old_RNA_TPM = old_data_TPM[, grep('RNA', colnames(old_data_TPM))]
colnames(old_Ribo_TPM)= paste0(substr(colnames(old_Ribo_TPM),1,5), isolates, '_0')
colnames(old_RNA_TPM)= paste0(substr(colnames(old_RNA_TPM),1,4), isolates, '_0')


new_TPM = fread('canonical_mapping/csv/counts.TPM.csv',data.table = F)
rownames(new_TPM)= new_TPM$V1
new_TPM= new_TPM[,-c(1,2,35,36)]
all_Ribo_TPM= cbind(new_TPM[intersect(rownames(new_TPM),rownames(old_Ribo_TPM)),grep('Ribo',colnames(new_TPM))], old_Ribo_TPM[intersect(rownames(new_TPM),rownames(old_Ribo_TPM)),])
all_Ribo_TPM=  all_Ribo_TPM[,order(colnames(all_Ribo_TPM))]

x = t(all_Ribo_TPM)
x=as.data.frame(x)
x$strain = substr(rownames(x),6,8 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,frame=T,label = TRUE)


all_RNA_TPM= cbind(new_TPM[intersect(rownames(new_TPM),rownames(old_RNA_TPM)),grep('RNA',colnames(new_TPM))], old_RNA_TPM[intersect(rownames(new_TPM),rownames(old_RNA_TPM)),])
all_RNA_TPM=  all_RNA_TPM[,order(colnames(all_RNA_TPM))]

x = t(all_RNA_TPM)
x=as.data.frame(x)
x$strain = substr(rownames(x),5,7 )

autoplot(prcomp(x[, -ncol(x)]), data = x, col = 'strain' ,frame=T,label = TRUE)



x= cor(all_Ribo_TPM, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)

x= cor(all_RNA_TPM, method = 's')
x=melt(x)
ggplot(x,aes(Var1, Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = 'lightblue',high = 'red',mid = 'white', midpoint = median(x$value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = round(value,2)), size = 4)


#### anota2Seq for data ####

anota2seq_pheno_vec = substr(colnames(RNA_seq),5,7)
ads <- anota2seqDataSetFromMatrix(
  dataP = ribo_seq,
  dataT = RNA_seq,
  phenoVec = anota2seq_pheno_vec,
  dataType = "RNAseq",
  normalize = T)

ads <- anota2seqRun(ads,onlyGroup = T)

#anota2seqPlotFC(ads, plotToFile = FALSE,selContrast =7)


x= apply(ads@dataT,1, var)
boxplot(apply(ads@dataT,1, var), apply(ads@dataP,1, var),log = 'y')
x= dist(t(ads@dataT))
x=as.matrix(x)[lower.tri(as.matrix(x))]
boxplot(as.matrix(dist(t(ads@dataP)))[lower.tri(as.matrix(dist(t(ads@dataP))))],x)

boxplot(cor(ads@dataT)[lower.tri(cor(ads@dataT))],cor(ads@dataP)[lower.tri(cor(ads@dataP))])
plot(cor(ads@dataT)[lower.tri(cor(ads@dataT))],cor(ads@dataP)[lower.tri(cor(ads@dataP))])

norm_RNA= ads@dataT
norm_Ribo= ads@dataP
rownames(norm_Ribo)==rownames(norm_RNA)
rownames(norm_Ribo)= rownames(new_TPM)[as.numeric(rownames(norm_Ribo))]
rownames(norm_RNA)= rownames(new_TPM)[as.numeric(rownames(norm_RNA))]


mean_norm_Ribo=lapply(isolates,function(i){
  i<<-i
  temp= norm_Ribo[,grep(i,colnames(norm_Ribo))]
  rowMeans(temp)
})
mean_norm_Ribo= do.call(cbind,mean_norm_Ribo)
colnames(mean_norm_Ribo)=isolates

mean_norm_RNA=lapply(isolates,function(i){
  i<<-i
  temp= norm_RNA[,grep(i,colnames(norm_RNA))]
  rowMeans(temp)
})
mean_norm_RNA= do.call(cbind,mean_norm_RNA)
colnames(mean_norm_RNA)=isolates


####DEG one vs all ####
df_strains = tibble(
  strain = isolates
  
  
)

pair_strains = tidyr::expand_grid(S1 = df_strains$strain, S2 = df_strains$strain) |>
  dplyr::mutate(pair = paste0(S1,'-',S2),
                is_identical = (S1 == S2),
                is_duplicated = find_duplicate_pairs(S1, S2))
sample_names = colnames(norm_Ribo)

df_group = tibble(sample = sample_names) |>
  separate(col  = sample, sep = '_', remove = F,
           into = c('type','strain','techrep')) |>
  left_join(df_strains,by='strain') |>
  mutate( strain  = factor(strain,levels = df_strains$strain),
          
  ) 
df_group$strain

#####  design matrix for experiment #####
design <- model.matrix(~0+df_group$strain)
colnames(design) <- df_strains$strain
print(design[,],max =50)

#####  fit gene linear model given design matrix Ribo #####
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
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)$estimate
fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate
x= data.frame(Odd_Ratio=c(fisher.test(e)$estimate,
                          fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate)
)
x$type= c('Metabolism','Other')
ggplot(x, aes(type,-1+Odd_Ratio, fill = type))+
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


ggplot(x[!x$Var1=='Total',], aes(fill = Var2, y=per, x= Var1))+
  geom_bar(stat = 'identity')+
  coord_polar(theta = "y", start = 0, direction = 1, clip = "on")+
  scale_fill_manual(values = c('gray23','gray71'))+
  theme_void()+
  geom_bar(data= x[x$Var1=='Total',], aes(fill = Var2, y=per, x= 2.6),stat = 'identity', width=0.2)

gostres <- gost(query = unique(one_vs_all_Ribo[one_vs_all_Ribo$is_DE,"gene"]) , organism = "scerevisiae", 
                ordered_query = FALSE,  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, user_threshold = 0.05, correction_method = "g_SCS",  
                domain_scope = "custom", custom_bg =rownames( norm_Ribo), numeric_ns = "", sources = NULL, as_short_link = FALSE) 
#view(gostres$result)


a = sum(unique(DE_gene)%in%comp)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_Ribo)%in%comp)-a
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)


a = sum(unique(DE_gene)%in%ess)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_Ribo)%in%ess)-a
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)



#####  fit gene linear model given design matrix RNA #####
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
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)$estimate
fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate
x= data.frame(Odd_Ratio=c(fisher.test(e)$estimate,
                          fisher.test(matrix(c(c,d,a,b),nrow= 2, byrow = T))$estimate)
)
x$type= c('Metabolism','Other')
ggplot(x, aes(type,-1+Odd_Ratio, fill = type))+
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


ggplot(x[!x$Var1=='Total',], aes(fill = Var2, y=per, x= Var1))+
  geom_bar(stat = 'identity')+
  coord_polar(theta = "y", start = 0, direction = 1, clip = "on")+
  scale_fill_manual(values = c('gray23','gray71'))+
  theme_void()+
  geom_bar(data= x[x$Var1=='Total',], aes(fill = Var2, y=per, x= 2.6),stat = 'identity', width=0.2)


gostres <- gost(query = unique(one_vs_all_RNA[one_vs_all_RNA$is_DE,"gene"]) , organism = "scerevisiae", 
                ordered_query = FALSE,  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, user_threshold = 0.05, correction_method = "g_SCS",  
                domain_scope = "custom", custom_bg =rownames( norm_RNA), numeric_ns = "", sources = NULL, as_short_link = FALSE) 
view(gostres$result)


a = sum(unique(DE_gene)%in%comp)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_Ribo)%in%comp)-a
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)


a = sum(unique(DE_gene)%in%ess)
b = length(unique(DE_gene))-a
c = sum(rownames(norm_Ribo)%in%ess)-a
d = length(rownames(prot_exp))-a-b-c
e = matrix(c(a,b,c,d),nrow= 2, byrow = T)
fisher.test(e)

#####real explo of Comp and ess in DEG#####
x = lapply(unique(one_vs_all_Ribo$isolate),function(i){
  temp= i
  i = one_vs_all_Ribo[one_vs_all_Ribo$isolate==i&one_vs_all_Ribo$is_DE,"gene"]
  
  a = sum(unique(i)%in%comp)
  b = length(unique(i))-a
  c = sum(rownames(norm_Ribo)%in%comp)-a
  d = length(rownames(prot_exp))-a-b-c
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
  d = length(rownames(prot_exp))-a-b-c
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
  
  a = sum(unique(i)%in%comp)
  b = length(unique(i))-a
  c = sum(rownames(norm_RNA)%in%comp)-a
  d = length(rownames(prot_exp))-a-b-c
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
  d = length(rownames(prot_exp))-a-b-c
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
deg_ess_comp[deg_ess_comp$p_val_corr>0.05,"data_plot"]=paste0('NA',deg_ess_comp[deg_ess_comp$p_val_corr>0.05,"data"])

ggplot(deg_ess_comp,aes(x= isolate, y =corrected, fill = data_plot))+
  geom_bar(stat = 'identity',position = 'dodge',col='black',alpha=1)+
  theme_classic()+
  facet_grid(~type)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5, 1), labels = c(0,0.5,1,1.5,2))+
  ylab('Odds ratios')+
  xlab('Isolates')+
  scale_fill_manual(values = c('lightgrey','lightgrey','lightsalmon','royalblue'))

#### Post-trans Buff exploration####

#####abs LFC#####
abs_log_RNA = pbapply(combn(isolates,2),2,function(i){ #compute abs lFC in each pairwise comparison for RNA
  temp= mean_norm_RNA[,c(grep(i[1], colnames(mean_norm_RNA)),grep(i[2], colnames(mean_norm_RNA)))]
  abs(log2(temp[,1]/temp[,2]))
  
})


abs_log_Ribo = pbapply(combn(isolates,2),2,function(i){#compute abs lFC in each pairwise comparison for Ribo
  temp= mean_norm_Ribo[,c(grep(i[1], colnames(mean_norm_Ribo)),grep(i[2], colnames(mean_norm_Ribo)))]
  abs(log2(temp[,1]/temp[,2]))
  
})

boxplot(melt(abs_log_RNA)[,3],melt(abs_log_Ribo)[,3],log='y')
wilcox.test(melt(abs_log_RNA)[,3],melt(abs_log_Ribo)[,3],log='y', paired = T)
mean(na.omit(melt(abs_log_RNA)[,3]))
mean(na.omit(melt(abs_log_Ribo)[,3]))

abs_lfc_delta = abs_log_Ribo-abs_log_RNA #difference between abs LFC Ribo and RNA (Ribo-RNA) if PTB, should be overall negative
median(na.omit(melt(abs_lfc_delta)[,3]))
temp = fgsea(pathways = GO2geneID,stats= sort(rowMeans(abs_lfc_delta,na.rm = T))) #any gene category more affected? --> GO
temp = temp[temp$padj<0.05,]
temp =temp[order(temp$pathway),]
temp$term= go2term(temp$pathway)[,2]

p.adjust(apply(abs_lfc_delta,1,function(i)(shapiro.test(na.omit(i))$p.value)),method = 'fdr') #are the 28 abs LFC value per gene normally distribuited

temp = rowMeans(abs_lfc_delta,na.rm = T)[p.adjust(apply(abs_lfc_delta,1,function(i)(wilcox.test(na.omit(i),mu = 0)$p.value)),method = 'fdr')<0.05] #since 28 are normally distribuited, t.test to check and select the confidently buffered or unbuffered genes 
median(temp)
temp=temp[temp<1]
names(temp)
buffered_genes = temp

#####dist #####

x= dist(t(mean_norm_Ribo))
x=as.matrix(x)[lower.tri(as.matrix(x))]
boxplot(as.matrix(dist(t(mean_norm_RNA)))[lower.tri(as.matrix(dist(t(mean_norm_RNA))))],x)
wilcox.test(as.matrix(dist(t(mean_norm_RNA)))[lower.tri(as.matrix(dist(t(mean_norm_RNA))))],x,paired = T)
median(as.matrix(dist(t(mean_norm_RNA)))[lower.tri(as.matrix(dist(t(mean_norm_RNA))))])
boxplot(as.matrix(dist(t(norm_RNA)))[lower.tri(as.matrix(dist(t(norm_RNA))))],as.matrix(dist(t(norm_Ribo)))[lower.tri(as.matrix(dist(t(norm_Ribo))))])
wilcox.test(as.matrix(dist(t(norm_RNA)))[lower.tri(as.matrix(dist(t(norm_RNA))))],as.matrix(dist(t(norm_Ribo)))[lower.tri(as.matrix(dist(t(norm_Ribo))))],paired=T)


#####correlation#####

x= cor(mean_norm_RNA, method = 's')[lower.tri(cor(mean_norm_RNA, method = 's'))]
temp= cor(mean_norm_Ribo, method = 's')[lower.tri(cor(mean_norm_Ribo, method = 's'))]
boxplot( x, temp)
wilcox.test(x, temp, paired = T)
median(x)
median(temp)

#####variance#####

temp = apply(mean_norm_Ribo,1,var)
x = apply(mean_norm_RNA,1,var)
boxplot(x,temp, log='y')
wilcox.test(x,temp, paired = T)
median(x)
median(temp)


#### anota2Seq exploration ####



#####generation of the pairwise comparison####
#annota_buff_each_strain= pblapply(isolates,function(i){
#  i<<-i
# 
#  anota2seq_pheno_vec = substr(colnames(RNA_seq),5,7)
#  anota2seq_pheno_vec[anota2seq_pheno_vec==i]=paste0(1,anota2seq_pheno_vec[anota2seq_pheno_vec==i])
#  ads <- anota2seqDataSetFromMatrix(
#   dataP = ribo_seq,
#   dataT = RNA_seq,
#   phenoVec = anota2seq_pheno_vec,
#   dataType = "RNAseq",
#   normalize = T)
# 
#  ads <- (anota2seqRun(ads,onlyGroup = T,useProgBar = F))
#  return(ads)
#})
#names(annota_buff_each_strain)=isolates
#write_rds(annota_buff_each_strain, 'buffering_anota2seq.RDS')
annota_buff_each_strain= readRDS('buffering_anota2seq.RDS')
annota_buff_each_strain$CQC@contrasts
annota_buff_each_strain$BPL@contrasts
anota2seqPlotFC(annota_buff_each_strain$CQC, plotToFile = FALSE,selContrast =4)
anota2seqPlotFC(annota_buff_each_strain$BPL, plotToFile = FALSE,selContrast =7)

##### what are the buffered genes####
######frequently buffered selection######

anota2seqPlotFC(annota_buff_each_strain$CQC, plotToFile = FALSE,selContrast =2)
anota2seqPlotFC(annota_buff_each_strain$BAN, plotToFile = FALSE,selContrast =7)
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

ggplot(temp,aes(ind, values))+
  geom_boxplot()+
  stat_compare_means()
range(temp$values)
mean(temp$values)

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
ggplot(temp,aes(strain, V1))+
  geom_boxplot()+
  stat_compare_means()

mean(unique(temp$V1))
range(temp$V1)

temp= stack(lapply(buffered_genes_anota, function(i){i<<-i
unlist(i)
}))
hist(table(temp$values)/2, xlab= 'Occurence in buffered genes', ylab ='Number of genes',main='')
unique(temp$values)
x = as.data.frame(table(temp$values)/2)
x$Var1= as.character(x$Var1)
x=rbind(x,data.frame(Var1=rownames(norm_Ribo)[!rownames(norm_Ribo)%in%x$Var1],Freq= rep(0,length(rownames(norm_Ribo)[!rownames(norm_Ribo)%in%x$Var1]))))
sum(x$Freq!=0)
sum(x$Freq>14)

x$freq_buff= x$Freq>14


######Ess or comp######
x$ess = x$Var1%in%ess
x$comp = x$Var1%in%comp
ggplot(x,aes(comp,Freq))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  stat_compare_means()

ggplot(x,aes(ess,Freq))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  stat_compare_means()



a = length(intersect(ess, x[x$Freq>14,"Var1"]))
b=  length(x[x$Freq>14,"Var1"])-a
c= length(intersect(ess, x[x$Freq<15,"Var1"]))
d= nrow(x)-a-b-c
fisher.test(matrix(c(a,b,c,d),nrow=2))

a = length(intersect(comp, x[x$Freq>14,"Var1"]))
b=  length(x[x$Freq>14,"Var1"])-a
c= length(intersect(comp, x[x$Freq<15,"Var1"]))
d= nrow(x)-a-b-c
fisher.test(matrix(c(a,b,c,d),nrow=2))

######functional enrichment #######

temp=setNames(x$Freq,nm = x$Var1)
temp=sort(temp)
temp= fgsea(temp,pathways=GO2geneID )

temp = na.omit(temp)
temp= temp[temp$padj<0.1,]
temp=temp[order(temp$pathway),]
temp$term=go2term(temp$pathway)[,2]
write_clip(temp[order(temp$NES),-c(1,2,4,5,8)])

######mRNA level######


plot(rowMeans(norm_RNA), apply(norm_RNA,1, var), log='y')
text(3,20,cor.test(rowMeans(norm_RNA), apply(norm_RNA,1, var), method = 's')$estimate)




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

ggplot(temp[!duplicated(temp$value),],aes(variable, value))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme_classic()+
  stat_compare_means()


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
ggplot(temp, aes(label, col = type, y = values))+
  geom_boxplot()+
  theme_classic()
#  stat_compare_means()

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

######new buffering of complex and ess####

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
#buff_comp[!buff_comp$sign,"OR_cor"]=NA
ggplot(buff_comp,aes(strain,vs, col=OR_cor,size= log10_pval))+
  geom_point()+
  theme_test()+
  scale_color_gradient2('Odds ratio',midpoint = 1,low = 'blue',high = 'red',mid = 'white',na.value = 'transparent',lim= c(0,4))+
  scale_size('-Log10(p-value)')+
  ylab('Isolates')+
  xlab('Isolates')





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
  scale_size('-Log10(p-value)')+
  ylab('Isolates')+
  xlab('Isolates')


######tAI index######


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
ggplot(temp,aes(variable, value))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme_classic()+
  stat_compare_means()

#### Exploration for Pangenome ####



rep_data_pan = fread("allORFs_pangenome/counts.allORFs_pangenome.csv", data.table = F, header = T)
rownames(rep_data_pan)=rep_data_pan$V1
RNA_seq_pan = rep_data_pan[,grep('RNA',colnames(rep_data_pan))]
ribo_seq_pan = rep_data_pan[,grep('Ribo',colnames(rep_data_pan))]
unique(substr(colnames(ribo_seq_pan),6,8))==unique(substr(colnames(RNA_seq_pan),5,7))
isolates = unique(substr(colnames(RNA_seq_pan),5,7))
anota2seq_pheno_vec = substr(colnames(RNA_seq_pan),5,7)
rownames(RNA_seq_pan)= rep_data_pan$V1
rownames(ribo_seq_pan)= rep_data_pan$V1

ads <- anota2seqDataSetFromMatrix(
  dataP = ribo_seq_pan,
  dataT = RNA_seq_pan,
  phenoVec = anota2seq_pheno_vec,
  dataType = "RNAseq",
  normalize = TRUE,filterZeroGenes =F)

norm_Ribo_pan= ads@dataP
norm_RNA_pan= ads@dataT
tpm_pangenome= cbind(norm_RNA_pan,norm_Ribo_pan)


#####number of expressed accessory#####

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


#####HGT BPL ######
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

temp= tpm_pangenome[,grep('BPL', colnames(tpm_pangenome))]
temp[rep_data_pan[rownames(temp),colnames(temp)]==0]=0
temp=temp[!rowSums(temp[,1:2])==0,]
temp=temp[!rowSums(temp[,3:4])==0,]
x=temp
temp=melt(temp)
temp$hgt=temp$Var1%in%rownames(b) 
sum(temp$hgt)


ggplot(temp,aes(x= Var2,y=value, fill= hgt))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c('grey', 'salmon'))

x=as.data.frame(x)
x$TE_1=x$Ribo_BPL_1/x$RNA_BPL_1

x$TE_2=x$Ribo_BPL_2/x$RNA_BPL_2
x$HGT = rownames(x)%in%rownames(b)


ggplot(x[!x$HGT,],aes(RNA_BPL_1,Ribo_BPL_1))+
  geom_point()+
  theme_classic()+
  geom_point(data=x[x$HGT,],aes(RNA_BPL_1,Ribo_BPL_1), col='magenta')+
  stat_cor(method = 's')+
  geom_abline(intercept =  0,slope =1)

ggplot(x[,c('TE_1','HGT')],aes(HGT, TE_1, fill=HGT))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means()+
  scale_fill_manual(values = c('grey', 'salmon'))

ggplot(x[,c('TE_2','HGT')],aes(HGT, TE_2, fill=HGT))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means()+
  scale_fill_manual(values = c('grey', 'salmon'))


#####HGT BPL reference normalized value######
temp= b[,grep('BPL',colnames(b))]
temp = rbind(temp,cbind(norm_RNA,norm_Ribo)[,grep('BPL',colnames(cbind(norm_RNA,norm_Ribo)))])
temp=as.data.frame(temp)
temp$HGT = rownames(temp)%in%rownames(b)
ggplot(melt(temp),aes(variable, value, col= HGT))+
  geom_boxplot()
temp$TE_1=temp$Ribo_BPL_1/temp$RNA_BPL_1

temp$TE_2=temp$Ribo_BPL_2/temp$RNA_BPL_2



ggplot(temp[!temp$HGT,],aes(RNA_BPL_1,Ribo_BPL_1))+
  geom_point()+
  theme_classic()+
  geom_point(data=temp[temp$HGT,],aes(RNA_BPL_1,Ribo_BPL_1), col='salmon')+
  stat_cor(method = 's')+
  geom_abline(intercept =  0,slope =1)

ggplot(temp[,c('TE_1','HGT')],aes(HGT, TE_1, fill=HGT))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means()+
  scale_y_log10()

ggplot(temp[,c('TE_2','HGT')],aes(HGT, TE_2, fill=HGT))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means()+
  scale_y_log10()


temp$TE_mean=rowMeans(temp[,c("TE_1","TE_2")])
ggplot(temp,aes(HGT, TE_mean, fill=HGT))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means()+
  scale_y_log10()

######tAI of these HGT vs the rest######
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
ggplot(temp,aes(Type, TE_mean, fill=Type))+
  geom_violin()+
  geom_boxplot(width=0.2,fill='white')+
  theme_classic()+
  stat_compare_means(label.x = 1.4,label = 'p.signif')+
  scale_y_log10()+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c('grey','#FAA094FF'))

table(temp$Type)
aggregate(data=temp, FUN = mean, TE_mean~Type)


x= rbind(data.frame(tAI=tai2_allORF_BPL[rownames(temp)[temp$Type=='HGT']], type='HGT'),data.frame(tAI=tai2_BPL, type='cannonical'))
ggplot(x,aes(type,tAI))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_classic()+
  stat_compare_means()

#####CPI introgression#####
het_orf = read.csv('Data/tpm_allInfo.csv')
all_var_gene <- readxl::read_xlsx('Data/8strainsSupGenes.xlsx')

all_var_gene <- as.data.frame(all_var_gene)
for( i in isolates){assign(paste('sup_ORF_',i,sep=''), all_var_gene[which(all_var_gene[,i]==1),])}


j = 5
b = sup_orf[,c(1,j+2)]
rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b
b = rep_data_pan[rep_data_pan$V1%in%x  ,grep('CPI',colnames(rep_data_pan))]
b = b[-which(rowSums(b)==0),]
b=b[!rowSums(b[,1:2])==0,]
b=b[!rowSums(b[,3:4])==0,]
b= tpm_pangenome[rownames(b),]
b[rep_data_pan[rownames(b),colnames(b)]==0]=0



c = rownames(b)[which(orf_info[rownames(b),7]=='Introgression')]
c = c[orf_info[c,3] %in% rownames(norm_Ribo)]
b= b[c,]

a = orf_info[c,3]
temp=c()
for(i in 1:length(a)){
  temp=c(temp, rownames(tpm_pangenome)[grep(a[i], rownames(tpm_pangenome))][1])
}
b = b[-which(is.na(temp)),]
x = het_orf[which(het_orf$Standardized_name=='CPI'),]  
x = x[which(x$systematic_name%in%a),]
x = x[which(x$AlleleVersion=='Sapa'),]

b= b[x$Annotation_Name,]
a= x$systematic_name
temp=c()
for(i in 1:length(a)){
  temp=c(temp, rownames(tpm_pangenome)[grep(a[i], rownames(tpm_pangenome))][1])
}
x = c()
for(i in isolates){
  if(i != 'CPI'){
    a = get(paste('sup_ORF_',i, sep=''))
    x=c(x,intersect(rownames(b),a$rn))
    print(i)}
}
if(length(x)>0){
  b = b[-which(rownames(b)==x),] 
  
  a = orf_info[x,3]
  a = rownames(tpm_pangenome)[grep(a, rownames(tpm_pangenome))]
  temp= temp[-which(is.na(temp))]}
x= tpm_pangenome[temp,]
x = x[,-grep('CPI',colnames(x))]  
x=cbind(x,b[,grep('CPI', colnames(b))])

introgression_data = x


######RNA######
a =introgression_data[,grep('RNA',colnames(introgression_data))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
p1=ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.25, paired = T)+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('log2TMM RNA-seq')+
  scale_fill_manual(values = c('grey','grey'))

######Ribo######

a =introgression_data[,grep('Ribo',colnames(introgression_data))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
p2 = ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.25, paired = T)+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('log2TMM Ribo-seq')+
  scale_fill_manual(values = c('grey','grey'))
grid.arrange(p1,p2,ncol=2)


#####CPI introgression with reference normalized values#####
het_orf = read.csv('Data/tpm_allInfo.csv')
all_var_gene <- readxl::read_xlsx('Data/8strainsSupGenes.xlsx')

all_var_gene <- as.data.frame(all_var_gene)
for( i in isolates){assign(paste('sup_ORF_',i,sep=''), all_var_gene[which(all_var_gene[,i]==1),])}


j = 5
b = sup_orf[,c(1,j+2)]
rownames(b)<- b$X
b <-rownames(b)[which(!is.na(b[,2]))]
x=b
b = rep_data_pan[rep_data_pan$V1%in%x  ,grep('CPI',colnames(rep_data_pan))]
b = b[-which(rowSums(b)==0),]
b=b[!rowSums(b[,1:2])==0,]
b=b[!rowSums(b[,3:4])==0,]
b= tpm_pangenome[rownames(b),]
b[rep_data_pan[rownames(b),colnames(b)]==0]=0



c = rownames(b)[which(orf_info[rownames(b),7]=='Introgression')]
c = c[orf_info[c,3] %in% rownames(norm_Ribo)]
b= b[c,]

a = orf_info[c,3]
temp=c()
for(i in 1:length(a)){
  temp=c(temp, rownames(tpm_pangenome)[grep(a[i], rownames(tpm_pangenome))][1])
}
b = b[-which(is.na(temp)),]
x = het_orf[which(het_orf$Standardized_name=='CPI'),]  
x = x[which(x$systematic_name%in%a),]
x = x[which(x$AlleleVersion=='Sapa'),]

b= b[x$Annotation_Name,]
a= x$systematic_name
temp=c()
for(i in 1:length(a)){
  temp=c(temp, rownames(tpm_pangenome)[grep(a[i], rownames(tpm_pangenome))][1])
}
x = c()
for(i in isolates){
  if(i != 'CPI'){
    a = get(paste('sup_ORF_',i, sep=''))
    x=c(x,intersect(rownames(b),a$rn))
    print(i)}
}
if(length(x)>0){
  b = b[-which(rownames(b)==x),] 
  
  a = orf_info[x,3]
  a = rownames(tpm_pangenome)[grep(a, rownames(tpm_pangenome))]
  temp= temp[-which(is.na(temp))]}
temp= substr(temp,6,12)
x= cbind(norm_RNA[temp,],norm_Ribo[temp,])
x = x[,-grep('CPI',colnames(x))]  
x=cbind(x,b[,grep('CPI', colnames(b))])


######RNA######
a =x[,grep('RNA',colnames(x))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
p1=ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.25, paired = T)+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('log2TMM RNA-seq')+
  scale_fill_manual(values = c('grey','grey'))

######Ribo######

a =x[,grep('Ribo',colnames(x))]

a = cbind(rowMeans(a[,-grep('CPI',colnames(a))]),rowMeans(a[,grep('CPI',colnames(a))]))
colnames(a)=c('Orthologs','Introgression')
a= melt(a)
p2 = ggplot(a, aes(x= Var2,y=value,fill= Var2))+
  geom_boxplot(width=0.3,outlier.colour = 'transparent')+
  geom_point(aes(fill = Var2), position = position_jitterdodge())+
  stat_compare_means(label.x = 1.25, paired = T)+
  theme_classic2()+
  theme(legend.position = 'none')+
  xlab('')+
  ylab('log2TMM Ribo-seq')+
  scale_fill_manual(values = c('grey','grey'))
grid.arrange(p1,p2,ncol=2)
