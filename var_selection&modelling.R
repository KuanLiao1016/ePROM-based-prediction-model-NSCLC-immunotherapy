library(mice)
library(survival)
library(dplyr)
library(rms)
library(StepReg)
library(lattice)
library(MASS)
library(SurvMetrics)


#select one imputed dataset and perform stepwise selection-----
imp_number <- 1:40
set.seed(315)
sample(imp_number,1) #imp number == 21

var_select_ds <- implong_pts_predict%>%filter(.imp == 21)
var_select_md0 <- coxph(Surv(surv_time,event)~ histology_coded + ps
                        + SmokingHistory + NLR + sex + firstline + PDL1Score + age,
                        var_select_ds)


summary(var_select_md0)
stepAIC(var_select_md0,direction = "backward")#ps, NLR, histology, age were selected, based on backward AIC
var_select_md1 <- coxph(Surv(surv_time/365.25, event) ~ histology_coded + ps + NLR  + rcs(age,3), var_select_ds)
var_select_md2 <- coxph(Surv(surv_time/365.25, event) ~ histology_coded + ps + NLR  + age, var_select_ds)
AIC(var_select_md1,var_select_md2)
var_select_md <- coxph(Surv(surv_time,event)~ histology_coded + ps
                        + NLR  + firstline + age,
                       var_select_ds)


#non-ePROM model-----
non_ePROM_md <- with(data=imp_pts_predict,
                     exp= coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + age))
summary(pool(non_ePROM_md))

res_md0 <- round(sapply(summary(pool(non_ePROM_md), conf.int = TRUE)[c(2,7,8)], exp),2)
res_md0 <- res_md0%>%as.data.frame()%>%mutate(summary(pool(non_ePROM_md))[c(1,6)],.before = 1)
res_md0

#eq_5d model----
eprom_md1 <- with(data=imp_pts_predict,
                  exp= coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + age + eq_5d))

res_md1 <- round(sapply(summary(pool(eprom_md1), conf.int = TRUE)[c(2,7,8)], exp),2)
res_md1 <- res_md1%>%as.data.frame()%>%mutate(summary(pool(eprom_md1))[c(1,6)],.before = 1)
res_md1
#symptom burden model----
eprom_md2 <- with(data=imp_pts_predict,
                  exp= coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + age+ symburden))

res_md2 <- round(sapply(summary(pool(eprom_md2), conf.int = TRUE)[c(2,7,8)], exp),2)
res_md2 <- res_md2%>%as.data.frame()%>%mutate(summary(pool(eprom_md2))[c(1,6)],.before = 1)
res_md2
#Moderate-to-severe symptom model----
eprom_md3 <- with(data=imp_pts_predict,
                  exp= coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + age + symcount))

res_md3 <- round(sapply(summary(pool(eprom_md3), conf.int = TRUE)[c(2,7,8)], exp),2)
res_md3 <- res_md3%>%as.data.frame()%>%mutate(summary(pool(eprom_md3))[c(1,6)],.before = 1)
res_md3




#lasso
library(glmnet)
tmp_predict <- var_select_ds%>%dplyr::select(SmokingHistory,PDL1Score,histology_coded,ps,age,NLR)
X <- model.matrix(~ . - 1, data = tmp_predict)
Y <- with(var_select_ds,Surv(surv_time, event))
lasso_md <- cv.glmnet(x=X,y = Y,alpha = 1, family = "cox", type.measure = "C")
plot(lasso_md)
lasso_md$lambda.min

