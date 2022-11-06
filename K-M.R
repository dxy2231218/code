rm(list = ls())
dir.create('Survival analysis')
library(pacman)
p_load(survival,survminer,dplyr,finalfit)
library(magrittr)

tcga_train <- read.table('model/risk.TCGAtrain.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga_train %<>% mutate(Risk=factor(Risk,levels = c('low','high'))) 
tcga_test <-read.table('model/risk.TCGAtest.txt',sep='\t',header = T,check.names = F,row.names = 1) 
tcga_test %<>% mutate(Risk=factor(Risk,levels = c('low','high')))
tcga_all<-read.table('model/risk.TCGAall.txt',sep='\t',header = T,check.names =F,row.names = 1)
tcga_all %<>% mutate(Risk=factor(Risk,levels = c('low','high')))
#geo<-read.table('model/risk.GEO.txt',sep='\t',header = T,check.names = F)
#geo %<>% mutate(Risk=factor(Risk,levels = c('low','high')))

do_surve<-function(dat){
  data<-get(dat)
  formula <- as.formula("Surv(futime,fustat)~Risk")
  fit <- surv_fit(formula,data = data)
  pval<-surv_pvalue(fit)[,2]
  if(pval<0.01){
    pval<-formatC(pval,digits = 2,format = 'E')
  }else pval<- as.character(signif(pval,3))
 
  dependent <- "Surv(futime,fustat==1)"
  explanory <- 'Risk'
  cox <- coxphuni(dependent = dependent,explanatory = explanory,.data = data)
  HR<-fit2df(cox,condense=F)[,2]
  HR<-signif(HR,3)
  ggsurvplot(fit)
  p = ggsurvplot(fit, data = data, xlab = "Time(Years)", ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                 size = 1.2, palette = c("navyblue", "firebrick1"), 
                 ylim = c(0, 101),
                 legend.title = 'Risk', legend.labs = c("low risk","high risk"), font.legend = 14, legend = c(0.8,0.8), 
                 pval = F, 
                 break.time.by = 2, 
                 risk.table = T, risk.table.height = 0.2, risk.table.fontsize = 5, risk.table.col = "strata", tables.y.text = F, 
                 font.x = 16, font.y = 16, font.xtickslab = 13, font.ytickslab = 13,
                 fun = function(y) y*100)  
                 ~italic(pvalue)~.(pval)),size=4)
  p$plot<-p$plot+annotate('text',x=3,y=5,label=bquote(HR~value~.(HR)),size=4)
  p$table<-p$table+ theme(axis.ticks = element_blank(), axis.line = element_blank(), 
                          axis.title.x = element_blank(), axis.text.x = element_blank(), 
                          axis.title.y = element_blank())
  #pdf('survival/TCGA_train_risk.pdf',height = 6,width = 6)
  pdf(paste0('Survival analysis/',dat,'_risk.pdf'),height = 6,width = 6)
  print(p)
  dev.off()
}
do_surve('tcga_all')
