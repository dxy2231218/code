
dir.create('GSEA')
rm(list = ls())
library(pacman)
setwd('R/project2/')
options(connectionObserver = NULL) 
devtools::install_github('GuangchuangYu/clusterProfiler') 
BiocManager::install('org.Hs.eg.db')
BiocManager::install('clusterProfiler')
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
p_load(tidyverse)


load('clean_rna_exp.Rdata')
table(str_sub(colnames(exp),1,12)) %>% table()
exp %<>% rename_all(~str_sub(.,1,12))

risk <- read.table('model/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)


common <- intersect(rownames(risk),colnames(exp))
risk <- risk[common,]
exp_data<-exp[,common]
high_risk_sample <- risk %>%filter(Risk=='high') %>% rownames()
low_risk_sample <- risk %>%filter(Risk=='low') %>% rownames()
exp_data_high <- exp_data[,high_risk_sample]
exp_data_low <- exp_data[,low_risk_sample]
mean_h <- rowMeans(exp_data_high)
mean_h[mean_h<0.00001]<-0.00001
mean_l <- rowMeans(exp_data_low)
mean_l[mean_l<0.00001]<-0.00001
logFC <- log2(mean_h)-log2(mean_l)
logFC<-sort(logFC,decreasing = T)

genes<- names(logFC)

gmt <- read.gmt('c2.cp.kegg.v7.4.symbols.gmt')


kk <- GSEA(logFC,TERM2GENE = gmt,pvalueCutoff = 1,eps = 0)
kktab <- as.data.frame(kk) %>% filter(pvalue<0.05)
write.table(kktab,file="GSEA/GSEA.result.txt",sep="\t",quote=F,row.names = F)
save(kk,file = 'GSEA/GSEA.result.Rdata')
load('GSEA.result.Rdata')

kkup <- kktab %>% filter(NES>0)
show_term <-rownames(kkup)[1:5]
gseaplot2(kk, show_term, base_size=8, title="Enriched in high risk group")
ggsave('GSEA.highRisk.pdf', width=7, height=5.5)


kkdown <- kktab %>% filter(NES<0)
show_term <-rownames(kkdown)[1:5]
gseaplot2(kk, show_term, base_size=8, title="Enriched in low risk group")
ggsave('GSEA.lowRisk.pdf', width=7, height=5.5)

kktab_sort<-kktab %>% dplyr::arrange(desc(abs(NES)))
