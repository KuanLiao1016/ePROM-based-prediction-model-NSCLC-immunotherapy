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
    brier(boot_md0, times = 365.25, newdata = resampled_data)
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
    brier(boot_md0, times = 365.25, newdata = impds)
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


summary(survfit(app_md,newdata = var_select_ds), times = 365.25)
riskRegression::predictRisk(app_md,newdata = var_select_ds, times = 365.25)


######################
#Overall performance #
######################
#Bootstrap Resampling
time1 <- Sys.time()
B <- 500
m <- 40
result <- matrix(nrow = B, ncol = 8) 
result_scaled <- matrix(nrow = B, ncol = 8) 
result_change <- data.frame()
result_change_scaled <- data.frame()
result_bs_app <- data.frame()
result_sbs_app <- data.frame()
summary(result_bs_app)
result_bs_corrected <- data.frame() 
result_sbs_corrected <- data.frame() 
result_change_bs_corrected <- data.frame() 
result_change_sbs_corrected <- data.frame() 
result_app <- data.frame() 
set.seed(1429)
for(i in 1:m) {
  impds <- implong_pts_predict %>% filter(.imp == m)
  climatological_fit <- survfit(Surv(surv_time, event) ~ 1, data = impds)
  climatological_bs <- mean(summary(climatological_fit, times = 365.25)$surv^2)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep = TRUE)
    samp_index[order(samp_index)]
    
    # Bootstrapping data
    resampled_data <- impds[samp_index, ]
    
    # Boot model
    boot_md0 <- coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + age, data = resampled_data)
    boot_md0_eq5d <- update(boot_md0, . ~ . + eq_5d)
    boot_md0_symb <- update(boot_md0, . ~ . + symburden)
    boot_md0_symc <- update(boot_md0, . ~ . + symcount)
    
    # Brier scores for resampled data
    bs_resample_0 <- brier(boot_md0, times = 365.25, newdata = resampled_data)
    bs_resample_eq5d <- brier(boot_md0_eq5d, times = 365.25, newdata = resampled_data)
    bs_resample_symb <- brier(boot_md0_symb, times = 365.25, newdata = resampled_data)
    bs_resample_symc <- brier(boot_md0_symc, times = 365.25, newdata = resampled_data)
    
    # Brier scores for original data
    bs_new <- brier(boot_md0, times = 365.25, newdata = impds)
    bs_new_eq5d <- brier(boot_md0_eq5d, times = 365.25, newdata = impds)
    bs_new_symb <- brier(boot_md0_symb, times = 365.25, newdata = impds)
    bs_new_symc <- brier(boot_md0_symc, times = 365.25, newdata = impds)
    
    #change in brier
    change_resample_eq5d <- bs_resample_eq5d$brier - bs_resample_0$brier
    change_resample_symb <- bs_resample_symb$brier - bs_resample_0$brier
    change_resample_symc <- bs_resample_symc$brier - bs_resample_0$brier
    
    
    # Change in Brier scores for original data
    change_new_eq5d <- bs_new_eq5d$brier - bs_new$brier
    change_new_symb <- bs_new_symb$brier - bs_new$brier
    change_new_symc <- bs_new_symc$brier - bs_new$brier
    
    # Optimism
    result[j, 1] <- bs_resample_0$brier - bs_new$brier
    result[j, 2] <- bs_resample_eq5d$brier - bs_new_eq5d$brier 
    result[j, 3] <- bs_resample_symb$brier -  bs_new_symb$brier
    result[j, 4] <- bs_resample_symc$brier - bs_new_symc$brier
    result[j, 5] <- bs_resample_0$brier
    result[j, 6] <- bs_resample_eq5d$brier
    result[j, 7] <- bs_resample_symb$brier
    result[j, 8] <- bs_resample_symc$brier
    
    # Store optimism for change in BS
    result_change[j, 1] <- change_new_eq5d - change_resample_eq5d
    result_change[j, 2] <- change_new_symb - change_resample_symb
    result_change[j, 3] <- change_new_symc - change_resample_symc
    # Calculate BSmax
    climatological_fit_resampled <- survfit(Surv(surv_time, event) ~ 1, data = resampled_data)
    climatological_bs_resampled <- mean(summary(climatological_fit, times = 365.25)$surv^2)
    
    # Calculate scaled Brier scores in the validation set
    scaled_bs_new <- 1- (bs_new$brier / climatological_bs) 
    scaled_bs_new_eq5d <- 1- (bs_new_eq5d$brier / climatological_bs) 
    scaled_bs_new_symb <- 1-(bs_new_symb$brier / climatological_bs) 
    scaled_bs_new_symc <- 1-(bs_new_symc$brier / climatological_bs) 
    # Calculate scaled Brier scores in the development set    
    scaled_bs_resample_0 <- 1-(bs_resample_0$brier / climatological_bs_resampled) 
    scaled_bs_resample_eq5d <- 1- (bs_resample_eq5d$brier / climatological_bs_resampled)
    scaled_bs_resample_symb <- 1-(bs_resample_symb$brier / climatological_bs_resampled) 
    scaled_bs_resample_symc <- 1-(bs_resample_symc$brier / climatological_bs_resampled)
    
    # Optimism for scaled Brier scores
    result_scaled[j, 1] <- scaled_bs_resample_0 - scaled_bs_new 
    result_scaled[j, 2] <- scaled_bs_resample_eq5d - scaled_bs_new_eq5d
    result_scaled[j, 3] <- scaled_bs_resample_symb - scaled_bs_new_symb
    result_scaled[j, 4] <- scaled_bs_resample_symc - scaled_bs_new_symc
    result_scaled[j, 5] <- scaled_bs_resample_0 
    result_scaled[j, 6] <- scaled_bs_resample_eq5d 
    result_scaled[j, 7] <- scaled_bs_resample_symb
    result_scaled[j, 8] <- scaled_bs_resample_symc
    # Change in scaled Brier scores for original data
    change_new_scaled_eq5d <- scaled_bs_new_eq5d - scaled_bs_new
    change_new_scaled_symb <- scaled_bs_new_symb - scaled_bs_new
    change_new_scaled_symc <- scaled_bs_new_symc - scaled_bs_new
    
    # Change in scaled Brier scores for resampled data
    change_resample_scaled_eq5d <- scaled_bs_resample_eq5d - scaled_bs_resample_0
    change_resample_scaled_symb <- scaled_bs_resample_symb - scaled_bs_resample_0
    change_resample_scaled_symc <- scaled_bs_resample_symc - scaled_bs_resample_0
    
    # Optimism for change in scaled Brier scores
    result_change_scaled[j, 1] <- change_new_scaled_eq5d - change_resample_scaled_eq5d
    result_change_scaled[j, 2] <- change_new_scaled_symb - change_resample_scaled_symb
    result_change_scaled[j, 3]  <- change_new_scaled_symc - change_resample_scaled_symc
  }
  
  boot_md1 <- coxph(Surv(surv_time, event) ~ histology_coded + ps + NLR + firstline + age, data = impds)
  boot_md1_eq5d <- update(boot_md1, . ~ . + eq_5d)
  boot_md1_symb <- update(boot_md1, . ~ . + symburden)
  boot_md1_symc <- update(boot_md1, . ~ . + symcount)

  
  bs_app <- brier(boot_md1, times = 365.25, newdata = impds)
  bs_app_eq5d <- brier(boot_md1_eq5d, times = 365.25, newdata = impds)
  bs_app_symb <- brier(boot_md1_symb, times = 365.25, newdata = impds)
  bs_app_symc <- brier(boot_md1_symc, times = 365.25, newdata = impds)
  
  
  # Calculate apparent scaled Brier scores 
  scaled_bs_app <- 1- (bs_app$brier / climatological_bs) 
  scaled_bs_app_eq5d <- 1- (bs_app_eq5d$brier / climatological_bs) 
  scaled_bs_app_symb <- 1- (bs_app_symb$brier / climatological_bs)
  scaled_bs_app_symc <- 1- (bs_app_symc$brier / climatological_bs)
  
  #Apparent bs
  result_bs_app[i,1] <- bs_app$brier
  result_bs_app[i,2:3] <- quantile(result[, 5], c(0.025, 0.975))
  result_bs_app[i,4] <- bs_app_eq5d$brier
  result_bs_app[i,5:6] <- quantile(result[, 6], c(0.025, 0.975))
  result_bs_app[i,7] <- bs_app_symb$brier
  result_bs_app[i,8:9] <- quantile(result[, 7], c(0.025, 0.975))
  result_bs_app[i,10] <- bs_app_symc$brier
  result_bs_app[i,11:12] <- quantile(result[, 8], c(0.025, 0.975))  
  
  #Apparent sbs
  result_sbs_app[i,1] <- scaled_bs_app
  result_sbs_app[i,2:3] <- quantile(result_scaled[, 5], c(0.025, 0.975))
  result_sbs_app[i,4] <- scaled_bs_app_eq5d
  result_sbs_app[i,5:6] <- quantile(result_scaled[, 6], c(0.025, 0.975))
  result_sbs_app[i,7] <- scaled_bs_app_symb
  result_sbs_app[i,8:9] <- quantile(result_scaled[, 7], c(0.025, 0.975))
  result_sbs_app[i,10] <- scaled_bs_app_symb
  result_sbs_app[i,11:12] <- quantile(result_scaled[, 8], c(0.025, 0.975)) 

  # Optimism corrected bs
  result_bs_corrected[i, 1:3] <- bs_app$brier - quantile(result[, 1], c(0.025, 0.5, 0.975))
  result_bs_corrected[i, 4:6] <- bs_app_eq5d$brier - quantile(result[, 2], c(0.025, 0.5, 0.975))
  result_bs_corrected[i, 7:9] <- bs_app_symb$brier - quantile(result[, 3], c(0.025, 0.5, 0.975))
  result_bs_corrected[i, 10:12] <- bs_app_symc$brier - quantile(result[, 4], c(0.025, 0.5, 0.975))
  
  # Optimism corrected sbs
  result_sbs_corrected[i, 1:3] <- scaled_bs_app - quantile(result_scaled[, 1], c(0.025, 0.5, 0.975))
  result_sbs_corrected[i, 4:6] <- scaled_bs_app_eq5d - quantile(result_scaled[, 2], c(0.025, 0.5, 0.975))
  result_sbs_corrected[i, 7:9] <- scaled_bs_app_symb - quantile(result_scaled[, 3], c(0.025, 0.5, 0.975))
  result_sbs_corrected[i, 10:12] <- scaled_bs_app_symc - quantile(result_scaled[, 4], c(0.025, 0.5, 0.975))
  

  # Calculate apparent change in Brier scores
  change_bs_app_eq5d <- bs_app_eq5d$brier - bs_app$brier
  change_bs_app_symb <- bs_app_symb$brier - bs_app$brier
  change_bs_app_symc <- bs_app_symc$brier - bs_app$brier
  
  # Optimism-corrected change in Brier scores
  result_change_bs_corrected[i, 1:3] <- change_bs_app_eq5d - quantile(result_change[, 1], c(0.025, 0.5, 0.975))
  result_change_bs_corrected[i, 4:6] <- change_bs_app_symb - quantile(result_change[, 2], c(0.025, 0.5, 0.975))
  result_change_bs_corrected[i, 7:9] <- change_bs_app_symc - quantile(result_change[, 3], c(0.025, 0.5, 0.975))
  
  # Calculate the apparent change in scaled Brier scores
  change_sbs_app_eq5d <- scaled_bs_app_eq5d - scaled_bs_app
  change_sbs_app_symb <- scaled_bs_app_symb - scaled_bs_app
  change_sbs_app_symc <- scaled_bs_app_symc - scaled_bs_app
  
  # Optimism-corrected change in scaled Brier scores
  result_change_sbs_corrected[i, 1:3] <- change_sbs_app_eq5d - quantile(result_change_scaled[, 1], c(0.025, 0.5, 0.975))
  result_change_sbs_corrected[i, 4:6] <- change_sbs_app_symb - quantile(result_change_scaled[, 2], c(0.025, 0.5, 0.975))
  result_change_sbs_corrected[i, 7:9] <- change_sbs_app_symc - quantile(result_change_scaled[, 3], c(0.025, 0.5, 0.975))
  
}
time2 <- Sys.time()
time2-time1
colnames(result_bs_app) <- c(
  "bs_app_mean", "bs_app_95CI_lower", "bs_app_95CI_upper", 
  "bs_app_eq5d_mean", "bs_app_eq5d_95CI_lower", "bs_app_eq5d_95CI_upper", 
  "bs_app_symb_mean", "bs_app_symb_95CI_lower", "bs_app_symb_95CI_upper", 
  "bs_app_symc_mean", "bs_app_95CI_lower", "bs_app_95CI_upper"
)

summary(result_bs_app)
colnames(result_sbs_app) <- c(
  "scaled_bs_app_mean", "scaled_bs_app_95CI_lower", "scaled_bs_app_95CI_upper", 
  "scaled_bs_app_eq5d_mean", "scaled_bs_app_eq5d_95CI_lower", "scaled_bs_app_eq5d_95CI_upper", 
  "scaled_bs_app_symb_mean", "scaled_bs_app_symb_95CI_lower", "scaled_bs_app_symb_95CI_upper", 
  "scaled_bs_app_symc_mean", "scaled_bs_app_symc_95CI_lower", "scaled_bs_app_symc_95CI_upper"
)
summary(result_sbs_app)

colnames(result_bs_corrected) <- c(
  "BS_Corrected_Lower", "BS_Corrected_Median", "BS_Corrected_Upper",
  "BS_Corrected_EQ5D_Lower", "BS_Corrected_EQ5D_Median", "BS_Corrected_EQ5D_Upper",
  "BS_Corrected_Symb_Lower", "BS_Corrected_Symb_Median", "BS_Corrected_Symb_Upper",
  "BS_Corrected_Symc_Lower", "BS_Corrected_Symc_Median", "BS_Corrected_Symc_Upper"
)
summary(result_bs_corrected)
colnames(result_sbs_corrected) <- c(
  "SBS_Corrected_Lower", "SBS_Corrected_Median", "SBS_Corrected_Upper",
  "SBS_Corrected_EQ5D_Lower", "SBS_Corrected_EQ5D_Median", "SBS_Corrected_EQ5D_Upper",
  "SBS_Corrected_Symb_Lower", "SBS_Corrected_Symb_Median", "SBS_Corrected_Symb_Upper",
  "SBS_Corrected_Symc_Lower", "SBS_Corrected_Symc_Median", "SBS_Corrected_Symc_Upper"
)
summary(result_sbs_corrected)
colnames(result_change_bs_corrected) <- c(
  "Change_BS_Corrected_EQ5D_Lower", "Change_BS_Corrected_EQ5D_Median", "Change_BS_Corrected_EQ5D_Upper",
  "Change_BS_Corrected_Symb_Lower", "Change_BS_Corrected_Symb_Median", "Change_BS_Corrected_Symb_Upper",
  "Change_BS_Corrected_Symc_Lower", "Change_BS_Corrected_Symc_Median", "Change_BS_Corrected_Symc_Upper"
)
summary(result_change_bs_corrected,digits = 2)
colnames(result_change_sbs_corrected) <- c(
  "Change_SBS_Corrected_EQ5D_Lower", "Change_SBS_Corrected_EQ5D_Median", "Change_SBS_Corrected_EQ5D_Upper",
  "Change_SBS_Corrected_Symb_Lower", "Change_SBS_Corrected_Symb_Median", "Change_SBS_Corrected_Symb_Upper",
  "Change_SBS_Corrected_Symc_Lower", "Change_SBS_Corrected_Symc_Median", "Change_SBS_Corrected_Symc_Upper"
)
summary(result_change_bs_corrected,digits = 2)


