library(ggrisk)
rm(list = ls())
dir.create('ggrisk')


train <- read.table('model/risk.TCGAtrain.txt',sep='\t',header = T,check.names = F,row.names = 1)
test <- read.table('model/risk.TCGAtest.txt',sep='\t',header = T,check.names = F,row.names = 1)
all <-read.table('model/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)

data_all<-all %>% select(-c(8,9))
data_train <- train %>% select(-c(8,9))
data_test<- test %>% select(-c(8,9))
fit_train <- cph(formula = Surv(futime,fustat==1)~DDN.AS1+DLEU1+RGS5+RUSC1.AS1+TMPO.AS1,data = data_train)
fit_test <- cph(formula = Surv(futime,fustat==1)~DDN.AS1+DLEU1+RGS5+RUSC1.AS1+TMPO.AS1,data = data_test)
fit_all <- cph(formula = Surv(futime,fustat==1)~DDN.AS1+DLEU1+RGS5+RUSC1.AS1+TMPO.AS1,data = data_all)
fit_train2<-coxph(formula = Surv(futime,fustat)~DDN.AS1+DLEU1+RGS5+RUSC1.AS1+TMPO.AS1,data = data_train)

p1<-ggrisk(fit_train,cutoff.show = F)
p2<-ggrisk(fit_test,cutoff.show = F)
p3<-ggrisk(fit_all,cutoff.show = F)
ggrisk(fit_train2)
ggsave('ggrisk/ggrisk_train.pdf',plot = p1,height = 6,width = 8)
ggsave('ggrisk/ggrisk_test.pdf',plot = p2,height = 6,width = 8)
ggsave('ggrisk/ggrisk_all.pdf',plot = p3,height = 6,width = 8)
