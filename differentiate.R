library(pacman)
library(magrittr)
p_load(tidyverse)
dir.create('project1')
dir.create('differentiate')
setwd('project1/')

miRNA<- read.table('symbol-countJ}>]>-IDW*;;.txt',sep='\t',header = T,check.names = F)

lncRNA <- read.table('lncRNA_gene_list',sep='\t')

lncRNA_exp <- miRNA %>% filter(id %in% lncRNA$V1)
lncRNA_exp %<>% mutate(id=str_replace(id,'-','\\.'))
lncRNA_exp %<>% rename_all(~str_sub(.,1,16))

lncRNA_exp %<>% filter(rowMeans(lncRNA_exp[,-1])>0)
table(table(lncRNA_exp$id))
lncRNA_exp %<>% plyr::ddply('id',plyr::colwise(mean))
save(lncRNA_exp,file = 'clean_lncRNA_exp.Rdata')

library(edgeR)
row.names(lncRNA_exp)<-lncRNA_exp$id
lncRNA_exp %<>% select(-1)
sampleLabel <- factor(c(rep('N',3),rep('T',306)),levels = c('N','T'))


dge <- DGEList(counts = lncRNA_exp,group = sampleLabel)

keep = rowSums(cpm(dge) > 0.5) > 3 
dge_filter = dge[keep, keep.lib.sizes = FALSE]
dim(dge_filter$counts) # [1] 18672   419


dge_norm = calcNormFactors(dge_filter)  # D,HO method = "TMM"


logcpm = cpm(dge_norm, log = TRUE, prior.count = 1)



design = model.matrix(~sampleLabel)
y = estimateDisp(dge_norm, design, robust = TRUE)


fit = glmQLFit(y, design)
res = glmQLFTest(fit)
result = topTags(res, n = Inf)$table


write_csv(data.frame(Symbol = rownames(result), result, check.names = F), "differentiate/lncRNA_DEA_result_info.csv")
write.csv(logcpm,'lncRNA_logcpm.csv')


DGE <- result %>% mutate(index=case_when(logFC > 1 & FDR < 0.05 ~'up',
                                   logFC< -1 & FDR< 0.05 ~'down',
                                   TRUE ~'Not Sig'))
DGE %>% count(index)

DGE %>% ggplot(aes(logFC,-log10(FDR)))+geom_point(aes(color=index))+
  theme_bw()+scale_y_continuous(limits=c(0,25))+
  scale_color_manual(name='',values = c('red','grey50','blue'))
ggsave('differentiate/vol.pdf',height = 6,width = 8)
dge_list<-DGE %>% filter(index != 'Not Sig') %>% row.names(
  .
)
save(dge_list,file = 'differentiate/DiffGenes.Rdata')
dge_up <- DGE %>% filter(index=='up') %>% rownames(.)
dge_down <- DGE %>% filter(index=='down') %>% rownames(.)
diff <- c(dge_up,dge_down)

diff <- logcpm[diff,]
names(sampleLabel)<-colnames(lncRNA_exp)
type<- as.data.frame(sampleLabel)
pdf('differentiate/heatmap.pdf',height = 10,width = 14)
pheatmap::pheatmap(diff, 
                   annotation=type, 
                   color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
                   cluster_cols =F,
                   cluster_rows = F,
                   scale="row",
                   show_colnames = F,
                   show_rownames = T,
                   fontsize = 8,
                   fontsize_row=3,
                   fontsize_col=8)
dev.off()
load('differentiate/DiffGenes.Rdata')


lncRNA_exp <- read.csv('differentiate/lncRNA_logcpm.csv',header = T,check.names = F,row.names = 1)

lncRNA_diff <- lncRNA_exp[dge_list,]
lncRNA_diff %<>% select_if(str_detect(colnames(lncRNA_diff),'01A$'))

exp <- data.table::fread('symbol.txt',sep='\t',header = T,check.names = F)
exp %<>% as.data.frame(exp)
exp %<>% mutate(id=str_replace_all(id,'-','\\.'))
exp %<>% rename_all(~str_sub(.,1,16))

exp %<>% filter(rowMeans(exp[,-1])>0.5)
table(table(exp$id))
dim(exp)
exp <- plyr::ddply(exp,'id',plyr::colwise(mean))
rownames(exp)<- exp$id
exp %<>% select(-1)
exp %<>% select_if(str_detect(colnames(exp),'01A$'))
dim(exp)
save(exp,file = 'clean_rna_exp.Rdata')


genes <- read.table('gene.txt',sep='\t',header = T)[,1]
common <- intersect(colnames(exp),colnames(lncRNA_diff))
data <- cbind(t(exp[genes,])[common,],t(lncRNA_diff)[common,])
data %<>% as.data.frame() %>% select_if(~!all(is.na(.)))
dim(data)
data <- log2(data+0.01)
corr <- cor(data)
dim(corr)
corr_test<-corrplot::cor.mtest(corr)
corr_test_p <- corr_test$p
corr2 <- abs(corr[1:65,!colnames(corr)%in%genes]) > 0.4
dim(corr2)
class(corr2)
related_lncRNA <- corr2[,apply(corr2,2,function(m)sum(m)>=1)] %>% colnames()

apply(corr2,2,function(m)which(m==T))
corr_test_p['MEG8',c(94,169,346,400,447)]
which(genes=='FLT3LG')
lncRNA_diff['FLT3LG',]
exp['FLT3LG',]
save(corr,corr_test,corr2,related_lncRNA,file = 'differentiate/corr_relatred.Rdata')
