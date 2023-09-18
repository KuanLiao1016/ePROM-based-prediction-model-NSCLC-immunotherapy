########load packages############
library(survival)
library(mice)
library(SurvMetrics)
library(rms)
library(dplyr)
library(pec)
library(boot)
library(pROC)

#######Discrimination###########
#####Non-ePROM model-----

#Bootstrap Resampling
time1 <- Sys.time()
B <- 500
m <- 40
result <- matrix(nrow = B,ncol = 6)
result_imp <- data.frame()
result_app <- data.frame()
set.seed(1429)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    samp_index[order(samp_index)]
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                                  NLR + age, data = resampled_data)
    
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    AUC_boot <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                           resampled_data, 
                                           reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                                       resampled_data, 
                                       reverse = TRUE,
                                       timewt = "n/G2")$concordance
    
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    AUC_orig <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                           impds, 
                                           reverse = TRUE)$concordance
    result[j,1]
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                       impds, 
                                       reverse = TRUE,
                                       timewt = "n/G2")$concordance
    #optimism
    result[j,1] <- AUC_boot$auc - AUC_orig$auc
    result[j,2] <-  harrell_C_boot - harrell_C_orig
    result[j,3] <-   Uno_C_boot - Uno_C_orig
    result[j,4] <- AUC_orig$auc
    result[j,5] <- harrell_C_orig
    result[j,6] <- Uno_C_orig
  }
  boot_md1 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                      NLR + firstline + age, data = impds)
  harrell_C0 <- concordance(boot_md1)$concordance
  app_hc_upci <-  concordance(boot_md1)$concordance + 1.96*sqrt(concordance(boot_md1)$var)
  app_hc_loci <- concordance(boot_md1)$concordance - 1.96*sqrt(concordance(boot_md1)$var)
  lp_app <- predict(boot_md1)
  pred_prob <- predict(boot_md1,type = "risk")
  Uno_C0 <- concordance(Surv(surv_time, event) ~ lp_app, 
                        impds, 
                        reverse = TRUE,
                        timewt = "n/G2")$concordance
  app_uno_upci <- Uno_C0 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  app_uno_loci <- Uno_C0 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  #apparent
  auc0<- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  app_auc_ci <- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$ci
  result_app[i,1:3] <- c(app_auc_ci)
  result_app[i,4:6] <- c(harrell_C0, app_hc_loci, app_hc_upci)
  result_app[i,7:9] <- c(Uno_C0,app_uno_loci,app_uno_upci)
  #optimism corrected, 95% CIs are based on location-shifted bootstrap method by Noma et al (2021)
  result_imp[i,1:3] <- app_auc_ci-mean(result[,1])
  result_imp[i,4] <- harrell_C0 - mean(result[,2])
  result_imp[i,5] <- app_hc_loci - mean(result[,2])
  result_imp[i,6] <- app_hc_upci - mean(result[,2])
  result_imp[i,7] <- Uno_C0 - mean(result[,3])
  result_imp[i,8] <- app_uno_loci - mean(result[,3])
  result_imp[i,9] <- app_uno_upci - mean(result[,3])
}
time2 <- Sys.time()
time2-time1
colnames(result_imp) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp,digits = 2)
colnames(result_app) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_app,digits = 2)

####EQ5D model----
result_imp_1 <- data.frame()
result_app_1 <- data.frame()
set.seed(1429)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    samp_index[order(samp_index)]
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + age + eq_5d, data = resampled_data)
    
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    AUC_boot <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    AUC_orig <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    result[j,1]
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    #optimism
    result[j,1] <- AUC_boot$auc - AUC_orig$auc
    result[j,2] <-  harrell_C_boot - harrell_C_orig
    result[j,3] <-   Uno_C_boot - Uno_C_orig
    result[j,4] <- AUC_orig$auc
    result[j,5] <- harrell_C_orig
    result[j,6] <- Uno_C_orig
  }
  boot_md1 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                      NLR + firstline + age + eq_5d, data = impds)
  harrell_C0 <- concordance(boot_md1)$concordance
  app_hc_upci <-  concordance(boot_md1)$concordance + 1.96*sqrt(concordance(boot_md1)$var)
  app_hc_loci <- concordance(boot_md1)$concordance - 1.96*sqrt(concordance(boot_md1)$var)
  lp_app <- predict(boot_md1)
  pred_prob <- predict(boot_md1,type = "risk")
  Uno_C0 <- concordance(Surv(surv_time, event) ~ lp_app, 
                        impds, 
                        reverse = TRUE,
                        timewt = "n/G2")$concordance
  app_uno_upci <- Uno_C0 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  app_uno_loci <- Uno_C0 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  #apparent
  auc0<- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  app_auc_ci <- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$ci
  result_app_1[i,1:3] <- c(app_auc_ci)
  result_app_1[i,4:6] <- c(harrell_C0, app_hc_loci, app_hc_upci)
  result_app_1[i,7:9] <- c(Uno_C0,app_uno_loci,app_uno_upci)
  #optimism corrected, 95% CIs are based on location-shifted bootstrap method by Noma et al (2021)
  result_imp_1[i,1:3] <- app_auc_ci-mean(result[,1])
  result_imp_1[i,4] <- harrell_C0 - mean(result[,2])
  result_imp_1[i,5] <- app_hc_loci - mean(result[,2])
  result_imp_1[i,6] <- app_hc_upci - mean(result[,2])
  result_imp_1[i,7] <- Uno_C0 - mean(result[,3])
  result_imp_1[i,8] <- app_uno_loci - mean(result[,3])
  result_imp_1[i,9] <- app_uno_upci - mean(result[,3])
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_1) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_1, digits = 2)
colnames(result_app_1) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_app_1,digits = 2)



#####symptom burden model----

time1 <- Sys.time()
result_imp_2 <- data.frame()
result_app_2 <- data.frame()
set.seed(1429)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    samp_index[order(samp_index)]
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + age + symburden, data = resampled_data)
    
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    app_cstat_model <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    test_cstat_model <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    result[j,1]
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    #optimism
    result[j,1] <- app_cstat_model$auc - test_cstat_model$auc
    result[j,2] <-  harrell_C_boot - harrell_C_orig
    result[j,3] <-   Uno_C_boot - Uno_C_orig
  }
  boot_md1 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                      NLR + firstline + age + symburden, data = impds)
  harrell_C0 <- concordance(boot_md1)$concordance
  app_hc_upci <-  concordance(boot_md1)$concordance + 1.96*sqrt(concordance(boot_md1)$var)
  app_hc_loci <- concordance(boot_md1)$concordance - 1.96*sqrt(concordance(boot_md1)$var)
  lp_app <- predict(boot_md1)
  pred_prob <- predict(boot_md1,type = "risk")
  Uno_C0 <- concordance(Surv(surv_time, event) ~ lp_app, 
                        impds, 
                        reverse = TRUE,
                        timewt = "n/G2")$concordance
  app_uno_upci <- Uno_C0 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  app_uno_loci <- Uno_C0 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  auc0<- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  app_auc_ci <- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$ci
  result_app_2[i,1:3] <- c(app_auc_ci)
  result_app_2[i,4:6] <- c(harrell_C0, app_hc_loci, app_hc_upci)
  result_app_2[i,7:9] <- c(Uno_C0,app_uno_loci,app_uno_upci)
  #optimism corrected, 95% CIs are based on location-shifted bootstrap method by Noma et al (2021)
  result_imp_2[i,1:3] <- app_auc_ci-mean(result[,1])
  result_imp_2[i,4] <- harrell_C0 - mean(result[,2])
  result_imp_2[i,5] <- app_hc_loci - mean(result[,2])
  result_imp_2[i,6] <- app_hc_upci - mean(result[,2])
  result_imp_2[i,7] <- Uno_C0 - mean(result[,3])
  result_imp_2[i,8] <- app_uno_loci - mean(result[,3])
  result_imp_2[i,9] <- app_uno_upci - mean(result[,3])
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_2) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_2,digits = 2)
colnames(result_app_2) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_app_2,digits = 2)



##################################
#Moderate-to-severe symptom model#
##################################

time1 <- Sys.time()
result_imp_3 <- data.frame()
result_app_3 <- data.frame()
set.seed(1429)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    samp_index[order(samp_index)]
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + age + symcount, data = resampled_data)
    
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    app_cstat_model <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    test_cstat_model <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    result[j,1]
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    #optimism
    result[j,1] <- app_cstat_model$auc - test_cstat_model$auc
    result[j,2] <-  harrell_C_boot - harrell_C_orig
    result[j,3] <-   Uno_C_boot - Uno_C_orig
  }
  boot_md1 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                      NLR + firstline + age + symcount, data = impds)
  harrell_C0 <- concordance(boot_md1)$concordance
  app_hc_upci <-  concordance(boot_md1)$concordance + 1.96*sqrt(concordance(boot_md1)$var)
  app_hc_loci <- concordance(boot_md1)$concordance - 1.96*sqrt(concordance(boot_md1)$var)
  lp_app <- predict(boot_md1)
  pred_prob <- predict(boot_md1,type = "risk")
  Uno_C0 <- concordance(Surv(surv_time, event) ~ lp_app, 
                        impds, 
                        reverse = TRUE,
                        timewt = "n/G2")$concordance
  app_uno_upci <- Uno_C0 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                             impds, 
                                                             reverse = TRUE,
                                                             timewt = "n/G2")$var)
  app_uno_loci <- Uno_C0 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app, 
                                                 impds, 
                                                 reverse = TRUE,
                                                 timewt = "n/G2")$var)
  auc0<- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  app_auc_ci <- roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$ci
  result_app_3[i,1:3] <- c(app_auc_ci)
  result_app_3[i,4:6] <- c(harrell_C0, app_hc_loci, app_hc_upci)
  result_app_3[i,7:9] <- c(Uno_C0,app_uno_loci,app_uno_upci)
  #optimism corrected, 95% CIs are based on location-shifted bootstrap method by Noma et al (2021)
  result_imp_3[i,1:3] <- app_auc_ci-mean(result[,1])
  result_imp_3[i,4] <- harrell_C0 - mean(result[,2])
  result_imp_3[i,5] <- app_hc_loci - mean(result[,2])
  result_imp_3[i,6] <- app_hc_upci - mean(result[,2])
  result_imp_3[i,7] <- Uno_C0 - mean(result[,3])
  result_imp_3[i,8] <- app_uno_loci - mean(result[,3])
  result_imp_3[i,9] <- app_uno_upci - mean(result[,3])
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_3) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
colnames(result_app_3) <- c("AUC_lo","AUC","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_3,digits = 2)
summary(result_app_3,digits = 2)

###Change in discrimination measures----
summary(result_imp_1 - result_imp, digits = 2)
summary(result_app_1 - result_app, digits = 2)
summary(result_imp_2 - result_imp)
summary(result_app_2 - result_app)

result_imp_3 <- result_imp_3[, c(1, 3, 2, 4, 6, 5,7,9,8)]
summary(result_imp_3 - result_imp)
colnames(result_imp_3)
colnames(result_imp)
summary(result_app_3 - result_app)
#######Calibration###########

###Non-ePROM model----

#Internal validation corrected performance - calibration
time1 <- Sys.time()
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res <- data.frame()
intercept_res <- data.frame()
slope_res <- data.frame()
slope_fixed_res <- data.frame()
set.seed(666)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
for (j in 1:B) {
  samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
  #bootstrapping data
  resampled_data <- impds[samp_index,]
  #boot_model
  boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                     NLR + firstline + age, data = resampled_data,  x = T)
  
  #fixed time point - 1 year
  obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                         data = impds), 
                 times = horizon)
  
  pred <- riskRegression::predictRisk(boot_md0,
                                      newdata = impds,
                                      times = horizon)
  lp1 <- predict(boot_md0, newdata = impds, type = "lp")#linear predictors
  OE <- (1 - obj$surv) / mean(pred)
  OE
  l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  o_e_ratio <- c(OE,l_ci,u_ci)
  oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
  
  gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
  res_slope_fixed <- summary(gval)$coefficients
  slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
  
  #time range mean calibration
  p <- predict(boot_md0, newdata = impds, type = "expected")
  logbase <- log(p) - lp1
  
  pfit1 <- glm(event ~ offset(log(p)), 
               family = poisson, 
               data = impds,
               subset = (p > 0))
  
  pfit <- glm(event ~ lp1 + offset(logbase), 
              family = poisson, 
              data = impds,
              subset= (p > 0))
  
  res_oe <- summary(pfit1)$coefficients
  int_summary <- rbind(res_oe,int_summary)#mean calibration
  
  res_slope <- summary(pfit)$coefficients[2,]
  slope_summary <- rbind(res_slope,slope_summary)
}
OE_res[i,1] <- mean(oe_res[,1])
OE_res[i,2] <- mean(oe_res[,2])
OE_res[i,3] <- mean(oe_res[,3])
intercept_res[i,1] <- round(exp(mean(int_summary$Estimate)),2)
intercept_res[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
intercept_res[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
slope_res[i,1] <- round(mean(slope_summary[,1]),2)
slope_res[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
slope_res[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
slope_fixed_res[i,1] <- round(mean(slope_fixed_summary[,1]),2)
slope_fixed_res[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
slope_fixed_res[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res,digits= 2)#fixed time point mean calibration
summary(intercept_res,digits= 2)#time range mean calibration - standardised mortality ratio
summary(slope_res,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res,digits = 2)#fixed time weak calibration - calibration slope
time2 <- Sys.time()
time2-time1



#apparent performance
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res <- data.frame()
intercept_res_app <- data.frame()
slope_res_app <- data.frame()
slope_fixed_res_app <- data.frame()
OE_res_app <- data.frame()
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  app_md<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                   NLR + firstline + age, data = impds,  x = T)
    
    #fixed time point - 1 year
    obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                           data = impds), 
                   times = horizon)
    
    pred <- riskRegression::predictRisk(app_md,
                                        newdata = impds,
                                        times = horizon)
    lp1 <- predict(app_md, newdata = impds, type = "lp")#linear predictors
    OE <- (1 - obj$surv) / mean(pred)
    OE
    l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    o_e_ratio <- c(OE,l_ci,u_ci)
    oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
    
    gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
    res_slope_fixed <- summary(gval)$coefficients
    slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
    
    #time range mean calibration
    p <- predict(app_md, newdata = impds, type = "expected")
    logbase <- log(p) - lp1
    
    pfit1 <- glm(event ~ offset(log(p)), 
                 family = poisson, 
                 data = impds,
                 subset = (p > 0))
    
    pfit <- glm(event ~ lp1 + offset(logbase), 
                family = poisson, 
                data = impds,
                subset= (p > 0))
    
    res_oe <- summary(pfit1)$coefficients
    int_summary <- rbind(res_oe,int_summary)#mean calibration
    
    res_slope <- summary(pfit)$coefficients[2,]
    slope_summary <- rbind(res_slope,slope_summary)
  OE_res_app[i,1] <- mean(oe_res[,1])
  OE_res_app[i,2] <- mean(oe_res[,2])
  OE_res_app[i,3] <- mean(oe_res[,3])
  intercept_res_app[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_app[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_app[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_app[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_app[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_app[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_app[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_app[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_app[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_app, digits = 2)#fixed time point mean calibration
summary(intercept_res_app)#time range mean calibration - standardised mortality ratio
summary(slope_res_app,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_app,digits = 2)#fixed time point weak calibration - calibration slope





###EQ-5D model###

#Internal validation corrected performance - calibration
time1 <- Sys.time()
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res_eq5d <- data.frame()
intercept_res_eq5d <- data.frame()
slope_res_eq5d <- data.frame()
slope_fixed_res_eq5d <- data.frame()
set.seed(666)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + eq_5d, data = resampled_data,  x = T)
    
    #fixed time point - 1 year
    obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                           data = impds), 
                   times = horizon)
    
    pred <- riskRegression::predictRisk(boot_md0,
                                        newdata = impds,
                                        times = horizon)
    lp1 <- predict(boot_md0, newdata = impds, type = "lp")#linear predictors
    OE <- (1 - obj$surv) / mean(pred)
    OE
    l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    o_e_ratio <- c(OE,l_ci,u_ci)
    oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
    
    gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
    res_slope_fixed <- summary(gval)$coefficients
    slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
    
    #time range mean calibration
    p <- predict(boot_md0, newdata = impds, type = "expected")
    logbase <- log(p) - lp1
    
    pfit1 <- glm(event ~ offset(log(p)), 
                 family = poisson, 
                 data = impds,
                 subset = (p > 0))
    
    pfit <- glm(event ~ lp1 + offset(logbase), 
                family = poisson, 
                data = impds,
                subset= (p > 0))
    
    res_oe <- summary(pfit1)$coefficients
    int_summary <- rbind(res_oe,int_summary)#mean calibration
    
    res_slope <- summary(pfit)$coefficients[2,]
    slope_summary <- rbind(res_slope,slope_summary)
  }
  OE_res_eq5d[i,1] <- mean(oe_res[,1])
  OE_res_eq5d[i,2] <- mean(oe_res[,2])
  OE_res_eq5d[i,3] <- mean(oe_res[,3])
  intercept_res_eq5d[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_eq5d[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_eq5d[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_eq5d[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_eq5d[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_eq5d[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_eq5d[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_eq5d[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_eq5d[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_eq5d,digits= 2)#fixed time point mean calibration
summary(intercept_res_eq5d,digits= 2)#time range mean calibration - standardised mortality ratio
summary(slope_res_eq5d,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_eq5d,digits = 2)#fixed time weak calibration - calibration slope
time2 <- Sys.time()
time2-time1




#apparent performance
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res <- data.frame()
intercept_res_app_eq5d <- data.frame()
slope_res_app_eq5d <- data.frame()
slope_fixed_res_app_eq5d <- data.frame()
OE_res_app_eq5d <- data.frame()
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  app_md<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                   NLR + firstline + age +eq_5d, data = impds,  x = T)
  
  #fixed time point - 1 year
  obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                         data = impds), 
                 times = horizon)
  
  pred <- riskRegression::predictRisk(app_md,
                                      newdata = impds,
                                      times = horizon)
  lp1 <- predict(app_md, newdata = impds, type = "lp")#linear predictors
  OE <- (1 - obj$surv) / mean(pred)
  OE
  l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  o_e_ratio <- c(OE,l_ci,u_ci)
  oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
  
  gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
  res_slope_fixed <- summary(gval)$coefficients
  slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
  
  #time range mean calibration
  p <- predict(app_md, newdata = impds, type = "expected")
  logbase <- log(p) - lp1
  
  pfit1 <- glm(event ~ offset(log(p)), 
               family = poisson, 
               data = impds,
               subset = (p > 0))
  
  pfit <- glm(event ~ lp1 + offset(logbase), 
              family = poisson, 
              data = impds,
              subset= (p > 0))
  
  res_oe <- summary(pfit1)$coefficients
  int_summary <- rbind(res_oe,int_summary)#mean calibration
  
  res_slope <- summary(pfit)$coefficients[2,]
  slope_summary <- rbind(res_slope,slope_summary)
  OE_res_app_eq5d[i,1] <- mean(oe_res[,1])
  OE_res_app_eq5d[i,2] <- mean(oe_res[,2])
  OE_res_app_eq5d[i,3] <- mean(oe_res[,3])
  intercept_res_app_eq5d[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_app_eq5d[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_app_eq5d[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_app_eq5d[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_app_eq5d[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_app_eq5d[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_app_eq5d[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_app_eq5d[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_app_eq5d[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_app_eq5d,digits = 2)#fixed time point mean calibration
summary(intercept_res_app_eq5d,digits = 2)#time range mean calibration - standardised mortality ratio
summary(slope_res_app_eq5d,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_app_eq5d,digits = 2)#fixed time point weak calibration - calibration slope






######################
#Symptom burden model#
######################
time1 <- Sys.time()
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res_symb <- data.frame()
intercept_res_symb <- data.frame()
slope_res_symb <- data.frame()
slope_fixed_res_symb <- data.frame()
set.seed(666)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + symburden, data = resampled_data,  x = T)
    
    #fixed time point - 1 year
    obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                           data = impds), 
                   times = horizon)
    
    pred <- riskRegression::predictRisk(boot_md0,
                                        newdata = impds,
                                        times = horizon)
    lp1 <- predict(boot_md0, newdata = impds, type = "lp")#linear predictors
    OE <- (1 - obj$surv) / mean(pred)
    OE
    l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    o_e_ratio <- c(OE,l_ci,u_ci)
    oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
    
    gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
    res_slope_fixed <- summary(gval)$coefficients
    slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
    
    #time range mean calibration
    p <- predict(boot_md0, newdata = impds, type = "expected")
    logbase <- log(p) - lp1
    
    pfit1 <- glm(event ~ offset(log(p)), 
                 family = poisson, 
                 data = impds,
                 subset = (p > 0))
    
    pfit <- glm(event ~ lp1 + offset(logbase), 
                family = poisson, 
                data = impds,
                subset= (p > 0))
    
    res_oe <- summary(pfit1)$coefficients
    int_summary <- rbind(res_oe,int_summary)#mean calibration
    
    res_slope <- summary(pfit)$coefficients[2,]
    slope_summary <- rbind(res_slope,slope_summary)
  }
  OE_res_symb[i,1] <- mean(oe_res[,1])
  OE_res_symb[i,2] <- mean(oe_res[,2])
  OE_res_symb[i,3] <- mean(oe_res[,3])
  intercept_res_symb[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_symb[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_symb[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_symb[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_symb[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_symb[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_symb[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_symb[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_symb[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_symb,digits = 2)#fixed time point mean calibration
summary(intercept_res_symb,digits = 2)#time range mean calibration - standardised mortality ratio
summary(slope_res_symb,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_symb,digits = 2)#fixed time weak calibration - calibration slope
time2 <- Sys.time()
time2-time1


#apparent performance
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res <- data.frame()
intercept_res_app_symb <- data.frame()
slope_res_app_symb <- data.frame()
slope_fixed_res_app_symb <- data.frame()
OE_res_app_symb <- data.frame()
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  app_md<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                   NLR + firstline + age + symburden, data = impds,  x = T)
  
  #fixed time point - 1 year
  obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                         data = impds), 
                 times = horizon)
  
  pred <- riskRegression::predictRisk(app_md,
                                      newdata = impds,
                                      times = horizon)
  lp1 <- predict(app_md, newdata = impds, type = "lp")#linear predictors
  OE <- (1 - obj$surv) / mean(pred)
  OE
  l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  o_e_ratio <- c(OE,l_ci,u_ci)
  oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
  
  gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
  res_slope_fixed <- summary(gval)$coefficients
  slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
  
  #time range mean calibration
  p <- predict(app_md, newdata = impds, type = "expected")
  logbase <- log(p) - lp1
  
  pfit1 <- glm(event ~ offset(log(p)), 
               family = poisson, 
               data = impds,
               subset = (p > 0))
  
  pfit <- glm(event ~ lp1 + offset(logbase), 
              family = poisson, 
              data = impds,
              subset= (p > 0))
  
  res_oe <- summary(pfit1)$coefficients
  int_summary <- rbind(res_oe,int_summary)#mean calibration
  
  res_slope <- summary(pfit)$coefficients[2,]
  slope_summary <- rbind(res_slope,slope_summary)
  OE_res_app_symb[i,1] <- mean(oe_res[,1])
  OE_res_app_symb[i,2] <- mean(oe_res[,2])
  OE_res_app_symb[i,3] <- mean(oe_res[,3])
  intercept_res_app_symb[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_app_symb[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_app_symb[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_app_symb[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_app_symb[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_app_symb[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_app_symb[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_app_symb[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_app_symb[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_app_symb,digits = 2)#fixed time point mean calibration
summary(intercept_res_app_symb)#time range mean calibration - standardised mortality ratio
summary(slope_res_app_symb,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_app_symb,digits = 2)#fixed time point weak calibration - calibration slope



##################################
#Moderate_to_severe symptom model#
##################################
time1 <- Sys.time()
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res_symc <- data.frame()
intercept_res_symc <- data.frame()
slope_res_symc <- data.frame()
slope_fixed_res_symc <- data.frame()
set.seed(666)
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    #bootstrapping data
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + symcount, data = resampled_data,  x = T)
    
    #fixed time point - 1 year
    obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                           data = impds), 
                   times = horizon)
    
    pred <- riskRegression::predictRisk(boot_md0,
                                        newdata = impds,
                                        times = horizon)
    lp1 <- predict(boot_md0, newdata = impds, type = "lp")#linear predictors
    OE <- (1 - obj$surv) / mean(pred)
    OE
    l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
    o_e_ratio <- c(OE,l_ci,u_ci)
    oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
    
    gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
    res_slope_fixed <- summary(gval)$coefficients
    slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
    
    #time range mean calibration
    p <- predict(boot_md0, newdata = impds, type = "expected")
    logbase <- log(p) - lp1
    
    pfit1 <- glm(event ~ offset(log(p)), 
                 family = poisson, 
                 data = impds,
                 subset = (p > 0))
    
    pfit <- glm(event ~ lp1 + offset(logbase), 
                family = poisson, 
                data = impds,
                subset= (p > 0))
    
    res_oe <- summary(pfit1)$coefficients
    int_summary <- rbind(res_oe,int_summary)#mean calibration
    
    res_slope <- summary(pfit)$coefficients[2,]
    slope_summary <- rbind(res_slope,slope_summary)
  }
  OE_res_symc[i,1] <- mean(oe_res[,1])
  OE_res_symc[i,2] <- mean(oe_res[,2])
  OE_res_symc[i,3] <- mean(oe_res[,3])
  intercept_res_symc[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_symc[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_symc[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_symc[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_symc[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_symc[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_symc[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_symc[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_symc[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_symc,digits = 2)#fixed time point mean calibration
summary(intercept_res_symc)#time range mean calibration - standardised mortality ratio
summary(slope_res_symc,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_symc,digits = 2)#fixed time weak calibration - calibration slope
time2 <- Sys.time()
time2-time1


#apparent performance
int_summary <- data.frame()
horizon <- 365.25
alpha <- 0.05
oe_res <- data.frame()
slope_summary <- data.frame()
slope_fixed_summary <- data.frame()
OE_res <- data.frame()
intercept_res_app_symc <- data.frame()
slope_res_app_symc <- data.frame()
slope_fixed_res_app_symc <- data.frame()
OE_res_app_symc <- data.frame()
for(i in 1:m ){
  impds <- implong_pts_predict%>%filter(.imp == i)
  app_md<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                   NLR + rcs(PDL1Score, 3) + rcs(age, 3) +symcount, data = impds,  x = T)
  
  #fixed time point - 1 year
  obj <- summary(survfit(Surv(surv_time, event) ~ 1, 
                         data = impds), 
                 times = horizon)
  
  pred <- riskRegression::predictRisk(app_md,
                                      newdata = impds,
                                      times = horizon)
  lp1 <- predict(app_md, newdata = impds, type = "lp")#linear predictors
  OE <- (1 - obj$surv) / mean(pred)
  OE
  l_ci <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  u_ci <- OE * exp(qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
  o_e_ratio <- c(OE,l_ci,u_ci)
  oe_res <- rbind(o_e_ratio,oe_res)#mean calibration
  
  gval <- coxph(Surv(surv_time, event) ~ lp1, data = impds)
  res_slope_fixed <- summary(gval)$coefficients
  slope_fixed_summary <- rbind(res_slope_fixed,slope_fixed_summary)#weak calibration
  
  #time range mean calibration
  p <- predict(app_md, newdata = impds, type = "expected")
  logbase <- log(p) - lp1
  
  pfit1 <- glm(event ~ offset(log(p)), 
               family = poisson, 
               data = impds,
               subset = (p > 0))
  
  pfit <- glm(event ~ lp1 + offset(logbase), 
              family = poisson, 
              data = impds,
              subset= (p > 0))
  
  res_oe <- summary(pfit1)$coefficients
  int_summary <- rbind(res_oe,int_summary)#mean calibration
  
  res_slope <- summary(pfit)$coefficients[2,]
  slope_summary <- rbind(res_slope,slope_summary)
  OE_res_app_symc[i,1] <- mean(oe_res[,1])
  OE_res_app_symc[i,2] <- mean(oe_res[,2])
  OE_res_app_symc[i,3] <- mean(oe_res[,3])
  intercept_res_app_symc[i,1] <- round(exp(mean(int_summary$Estimate)),2)
  intercept_res_app_symc[i,2] <- round(exp(mean(int_summary$Estimate) + 1.96 * mean(int_summary$`Std. Error`)),2)
  intercept_res_app_symc[i,3] <- round(exp(mean(int_summary$Estimate) - 1.96 * mean(int_summary$`Std. Error`)),2)
  slope_res_app_symc[i,1] <- round(mean(slope_summary[,1]),2)
  slope_res_app_symc[i,2] <- round(mean(slope_summary[,1]) + 1.96 * mean(slope_summary[,2]),2)
  slope_res_app_symc[i,3] <- round(mean(slope_summary[,1]) - 1.96 * mean(slope_summary[,2]),2)
  slope_fixed_res_app_symc[i,1] <- round(mean(slope_fixed_summary[,1]),2)
  slope_fixed_res_app_symc[i,2] <- round(mean(slope_fixed_summary[,1]) + 1.96 * mean(slope_fixed_summary[,3]),2)
  slope_fixed_res_app_symc[i,3] <- round(mean(slope_fixed_summary[,1]) - 1.96 * mean(slope_fixed_summary[,3]),2)
}
summary(OE_res_app_symc,digits = 2)#fixed time point mean calibration
summary(intercept_res_app_symc)#time range mean calibration - standardised mortality ratio
summary(slope_res_app_symc,digits= 2)#time range weak calibration - calibration slope
summary(slope_fixed_res_app_symc,digits = 2)#fixed time point weak calibration - calibration slope







