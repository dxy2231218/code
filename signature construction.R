
dir.create('model')
rm(list = ls())
library(pacman)
library(magrittr)
p_load(tidyverse)
p_load(survival,survminer,glmnet,caret,timeROC,finalfit)

##load('clean_lncRNA_exp.Rdata')
lncRNA_exp<-read.table('lncRNA_logcpm.csv',sep=',',header = T,check.names = F)

related_lncRNA <- read.csv('gene_lncRNA correlation.csv',header = T)
related_lncRNA<-related_lncRNA %>% distinct(gene) %>% pull(gene)
colnames(lncRNA_exp)[1]<-'id'
lncRNA <- lncRNA_exp %>% filter(id %in% related_lncRNA)%>% column_to_rownames('id')


surv<-read.table('time.txt',sep='\t',check.names = F,header = T,row.names = 1)


lncRNA %<>% select_if(str_detect(colnames(lncRNA),'01A')) %>%
  rename_all(~str_sub(.,1,12))
common <- intersect(rownames(surv),colnames(lncRNA))
data <- cbind(surv[common,],t(lncRNA)[common,])
dim(data)
data %<>% mutate(futime=futime/365)
data %<>% rename_all(~str_replace_all(.,'-','\\.'))
do_model<-function(){
  intrain <- createDataPartition(y = data[,1],p = 0.65,list = F)
  train<-data[intrain,]
  test<- data[-intrain,]
  train_out <- train %>% mutate(id=rownames(train)) %>% relocate(id)
  test_out <- test %>% mutate(id=rownames(test)) %>% relocate(id)
  
 
  fit <- coxphuni(.data = train,dependent = "Surv(futime,fustat)",explanatory = setdiff(colnames(train),c('futime','fustat')))
  cox_res <- fit2df(fit,condense=F)

  out_uni_tab <- cox_res %>% filter(p<0.05)
  

  uni_sig_exp <- train_out %>% select(futime,fustat,out_uni_tab$explanatory)
  uni_sig_exp_out <- train_out %>% select(id,futime,fustat,out_uni_tab$explanatory)
  
 
  x<- as.matrix(uni_sig_exp[,3:ncol(uni_sig_exp)])
  y<-data.matrix(Surv(uni_sig_exp$futime,uni_sig_exp$fustat))
 
  fit <- glmnet(x,y,family = 'cox',maxit = 1000,alpha = 1) 
  cvfit <-cv.glmnet(x,y,family = 'cox',maxit = 1000) 
  #plot(fit)
  #plot(cvfit)
  coef <- coef(fit, s = cvfit$lambda.min) 
  index <- which(coef != 0)
  actCoef <- coef[index]
  lasso_gene <-rownames(coef)[index]
  lasso_sig_exp <- uni_sig_exp %>%select(futime,fustat,lasso_gene)
  lasso_sig_exp_out<- uni_sig_exp_out %>% select(id,futime,fustat,lasso_gene)
  gene_coef=cbind(Gene=lasso_gene, Coef=actCoef) 
  
  
  
  multicox_fit <- coxphmulti(.data = lasso_sig_exp,dependent = "Surv(futime,fustat)",explanatory = setdiff(colnames(lasso_sig_exp),c('futime','fustat')))
  multicox_fit2 <- coxph(Surv(futime, fustat) ~ ., data = lasso_sig_exp)
  multicox_fit2<- step(multicox_fit2,direction = 'both')
  multicox_res <- multicox_fit %>% fit2df(condense=F)
  multicox_sum <- summary(multicox_fit2)
  
  
  
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multicox_sum$coefficients[,"coef"],
    HR=multicox_sum$conf.int[,"exp(coef)"],
    HR.95L=multicox_sum$conf.int[,"lower .95"],
    HR.95H=multicox_sum$conf.int[,"upper .95"],
    pvalue=multicox_sum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  

  risk_score <- predict(multicox_fit2,newdata = train,type = 'risk') 
  cox_gene <- rownames(multicox_sum$coefficients)
  out_col <- c('futime','fustat',cox_gene)
  median_train_risk <- median(risk_score)
  risk <- ifelse(risk_score>median_train_risk,'high','low')
  train_risk_out <- cbind(id=rownames(train),train[,out_col],risk_score,Risk=risk)
  
  
  risk_score_test=predict(multicox_fit2,type="risk",newdata=test)
  risk_test <- ifelse(risk_score_test>median_train_risk,'high','low')
  test_risk_out <- cbind(id=rownames(test),test[,out_col],risk_score_test,Risk=risk_test)
  
  
  
 
  formula <- as.formula('Surv(futime,fustat==1)~Risk')
  fit_train <- surv_fit(formula,data =train_risk_out)
  ggsurvplot(fit_train)
  pvalue_train <- surv_pvalue(fit_train)[,2]
  fit_test <- surv_fit(formula,data = test_risk_out)
  #ggsurvplot(fit_test)
  pvalue_test <- surv_pvalue(fit_test)[,2]
  #fit_geo <- surv_fit(formula,data = geo_risk_out)
  #ggsurvplot(fit_geo)
  #pvalue_geo <- surv_pvalue(fit_geo)[,2]
  
  
  predict_time <- 1
  roc<-timeROC(T=train$futime, delta=train$fustat,
               marker=risk_score, cause=1,
               times=c(predict_time), ROC=TRUE)
  rocTest<-timeROC(T=test$futime, delta=test$fustat,
                   marker=risk_score_test, cause=1,
                   times=c(predict_time), ROC=TRUE)
  return(list(train_out=train_out,test_out=test_out,out_uni_tab=out_uni_tab,uni_sig_exp_out=uni_sig_exp_out,
              lasso_sig_exp_out=lasso_sig_exp_out,outMultiTab,train_risk_out,test_risk_out,
              fit_train,pvalue_train,fit_test,pvalue_test,roc,rocTest))
}
model1<-do_model()
while(T){
  model<-do_model()
  if((model[[10]])<0.05 && (model[[12]])<0.05){
    break
  }
}
model[[10]]
model[[12]]
model[[13]]$AUC[2]
model[[14]]$AUC[2]
model[[6]]
ggsurvplot(model[[9]])
model2<-model
save(model,file = 'model/model.Rdata')


write.table(model[[1]],file="modeldata.train.txt",sep="\t",quote=F,row.names=F)
write.table(model[[2]],file="model/data.test.txt",sep="\t",quote=F,row.names=F)


write.table(model[[3]],file="model/uni.trainCox.txt",sep="\t",row.names=F,quote=F)
write.table(model[[4]],file="model/uni.SigExp.txt",sep="\t",row.names=F,quote=F)


write.table(model[[5]],file="model/lasso.SigExp.txt",sep="\t",row.names=F,quote=F)


outMultiTab=model[[6]][,1:2]
write.table(outMultiTab,file="model/multiCox.txt",sep="\t",row.names=F,quote=F)
write.table(model[[7]],file="model/risk.TCGAtrain.txt",sep="\t",quote=F,row.names=F)
write.table(model[[8]],file="model/risk.TCGAtest.txt",sep="\t",quote=F,row.names=F)
#write.table(model[[9]],file="model/risk.GEO.txt",sep="\t",quote=F,row.names=F)



colnames(model[[7]])
colnames(model[[8]])
model[[8]] <- model[[8]] %>% rename(c('risk_score'='risk_score_test')) 
allRiskOut=rbind(model[[7]], model[[8]])
write.table(allRiskOut,file="model/risk.TCGAall.txt",sep="\t",quote=F,row.names=F)


uni_sig_exp <- read.table('model/uni.SigExp.txt',sep='\t',check.names = F,header = T,row.names = 1)
x<- as.matrix(uni_sig_exp[,3:ncol(uni_sig_exp)])
y<-data.matrix(Surv(uni_sig_exp$futime,uni_sig_exp$fustat))
fit <- glmnet(x,y,family = 'cox',maxit = 1000,alpha = 1)
cvfit <-cv.glmnet(x,y,family = 'cox',maxit = 1000) 
pdf("model/lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
pdf("model/lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
