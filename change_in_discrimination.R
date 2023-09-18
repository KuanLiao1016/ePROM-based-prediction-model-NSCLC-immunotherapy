library(survival)
library(SurvMetrics)
library(dplyr)
library(pec)
library(pROC)

####moderate-to-severe symptom model vs non-eprom model
time1 <- Sys.time()
B <- 500
m <- 40
result <- matrix(nrow = B,ncol = 6)
result_imp_c_3 <- data.frame()
result_app_c_3 <- data.frame()
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
    boot_md1 <- update(boot_md0, . ~ . + symcount)
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    lp_bs_1<- predict(boot_md1, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob_1 <- predict(boot_md1, type = "risk", time = 365.25)
    AUC_boot <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    AUC_boot_1 <- roc((event)~pred_prob_1,ci=TRUE, data = resampled_data)
    
    
    harrell_C_boot_1 <- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                                    resampled_data, 
                                    reverse = TRUE)$concordance
    
    Uno_C_boot_1<- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                               resampled_data, 
                               reverse = TRUE,
                               timewt = "n/G2")$concordance
    
    AUC_change <- AUC_boot_1$auc- AUC_boot$auc
    harrell_C_boot_change <- harrell_C_boot_1 -harrell_C_boot
    Uno_C_boot_change <- Uno_C_boot_1- Uno_C_boot
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    lp_orig_1 <- predict(boot_md1, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    pred_prob_vl_1 <- predict(boot_md1, type = "risk", time = 365.25, newdata = impds)
    AUC_orig <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    AUC_orig_1 <- roc((event)~pred_prob_vl_1,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    harrell_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                    impds, 
                                    reverse = TRUE)$concordance
    
    Uno_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                impds, 
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
    AUC_orig_change <- AUC_orig_1$auc - AUC_orig$auc
    harrell_C_orig_change <- harrell_C_orig_1- harrell_C_orig
    Uno_C_orig_change <- Uno_C_orig_1 - Uno_C_orig
    #optimism
    result[j,1] <- app_cstat_model_change - test_cstat_model_change
    result[j,2] <-  harrell_C_boot_change - harrell_C_orig_change
    result[j,3] <-   Uno_C_boot_change - Uno_C_orig_change
    result[j,4] <- AUC_orig_change
    result[j,5] <- harrell_C_orig_change
    result[j,6] <- Uno_C_orig_change
  }
  boot_md00 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                       NLR + firstline + age, data = impds)
  boot_md10 <- update(boot_md00, . ~ . + symcount)
  lp_app <- predict(boot_md00)
  lp_app_1 <- predict(boot_md10)
  pred_prob <- predict(boot_md00,type = "risk")
  pred_prob_1 <- predict(boot_md10,type = "risk")
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
  Uno_C0_1 <- concordance(Surv(surv_time, event) ~ lp_app_1, 
                          impds, 
                          reverse = TRUE,
                          timewt = "n/G2")$concordance
  
  app_uno_upci_1 <- Uno_C0_1 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  app_uno_loci_1 <- Uno_C0_1 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  
  Uno_C0_change <- Uno_C0_1 - Uno_C0
  harrell_C0_change <- concordance(boot_md10)$concordance -  concordance(boot_md00)$concordance
  auc0_change<- roc((surv_time<365.25)~pred_prob_1,ci=TRUE, data = impds)$auc - roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc

  #optimism corrected
  result_imp_c_3[i,1] <- auc0_change-mean(result[,1])
  result_imp_c_3[i,2:3] <- quantile(result[,4], c(0.025,0.975))
  result_imp_c_3[i,4] <- harrell_C0_change - mean(result[,2])
  result_imp_c_3[i,5:6] <- quantile(result[,5], c(0.025,0.975))
  result_imp_c_3[i,7] <- Uno_C0_change - mean(result[,3])
  result_imp_c_3[i,8:9] <- quantile(result[,6], c(0.025,0.975))
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_c_3) <- c("AUC","AUC_lo","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_c_3)


####symptom burden model vs non-eprom model
time1 <- Sys.time()
B <- 500
m <- 40
result <- matrix(nrow = B,ncol = 6)
result_imp_c_2 <- data.frame()
result_app_c_2 <- data.frame()
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
    boot_md1 <- update(boot_md0, . ~ . + symburden)
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    lp_bs_1<- predict(boot_md1, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob_1 <- predict(boot_md1, type = "risk", time = 365.25)
    AUC_boot <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    AUC_boot_1 <- roc((event)~pred_prob_1,ci=TRUE, data = resampled_data)
    
    
    harrell_C_boot_1 <- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                                    resampled_data, 
                                    reverse = TRUE)$concordance
    
    Uno_C_boot_1<- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                               resampled_data, 
                               reverse = TRUE,
                               timewt = "n/G2")$concordance
    
    AUC_change <- AUC_boot_1$auc- AUC_boot$auc
    harrell_C_boot_change <- harrell_C_boot_1 -harrell_C_boot
    Uno_C_boot_change <- Uno_C_boot_1- Uno_C_boot
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    lp_orig_1 <- predict(boot_md1, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    pred_prob_vl_1 <- predict(boot_md1, type = "risk", time = 365.25, newdata = impds)
    AUC_orig <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    AUC_orig_1 <- roc((event)~pred_prob_vl_1,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    harrell_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                    impds, 
                                    reverse = TRUE)$concordance
    
    Uno_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                impds, 
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
    AUC_orig_change <- AUC_orig_1$auc - AUC_orig$auc
    harrell_C_orig_change <- harrell_C_orig_1- harrell_C_orig
    Uno_C_orig_change <- Uno_C_orig_1 - Uno_C_orig
    #optimism
    result[j,1] <- app_cstat_model_change - test_cstat_model_change
    result[j,2] <-  harrell_C_boot_change - harrell_C_orig_change
    result[j,3] <-   Uno_C_boot_change - Uno_C_orig_change
    result[j,4] <- AUC_orig_change
    result[j,5] <- harrell_C_orig_change
    result[j,6] <- Uno_C_orig_change
  }
  boot_md00 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                       NLR + firstline + age, data = impds)
  boot_md10 <- update(boot_md00, . ~ . + symburden)
  lp_app <- predict(boot_md00)
  lp_app_1 <- predict(boot_md10)
  pred_prob <- predict(boot_md00,type = "risk")
  pred_prob_1 <- predict(boot_md10,type = "risk")
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
  Uno_C0_1 <- concordance(Surv(surv_time, event) ~ lp_app_1, 
                          impds, 
                          reverse = TRUE,
                          timewt = "n/G2")$concordance
  
  app_uno_upci_1 <- Uno_C0_1 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  app_uno_loci_1 <- Uno_C0_1 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  
  Uno_C0_change <- Uno_C0_1 - Uno_C0
  harrell_C0_change <- concordance(boot_md10)$concordance -  concordance(boot_md00)$concordance
  auc0_change<- roc((surv_time<365.25)~pred_prob_1,ci=TRUE, data = impds)$auc - roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  
  #optimism corrected
  result_imp_c_2[i,1] <- auc0_change-mean(result[,1])
  result_imp_c_2[i,2:3] <- quantile(result[,4], c(0.025,0.975))
  result_imp_c_2[i,4] <- harrell_C0_change - mean(result[,2])
  result_imp_c_2[i,5:6] <- quantile(result[,5], c(0.025,0.975))
  result_imp_c_2[i,7] <- Uno_C0_change - mean(result[,3])
  result_imp_c_2[i,8:9] <- quantile(result[,6], c(0.025,0.975))
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_c_2) <- c("AUC","AUC_lo","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_c_2)

####EQ-5D vs non-ePROM

time1 <- Sys.time()
B <- 500
m <- 40
result <- matrix(nrow = B,ncol = 6)
result_imp_c_1 <- data.frame()
result_app_c_1 <- data.frame()
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
    boot_md1 <- update(boot_md0, . ~ . + eq_5d)
    lp_bs<- predict(boot_md0, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob <- predict(boot_md0, type = "risk", time = 365.25)
    lp_bs_1<- predict(boot_md1, type = 'lp') # predict lp from the bootstrap model in the bs sample
    pred_prob_1 <- predict(boot_md1, type = "risk", time = 365.25)
    AUC_boot <- roc((event)~pred_prob,ci=TRUE, data = resampled_data)
    
    harrell_C_boot <- concordance(Surv(surv_time, event) ~ lp_bs, 
                                  resampled_data, 
                                  reverse = TRUE)$concordance
    
    Uno_C_boot<- concordance(Surv(surv_time, event) ~ lp_bs, 
                             resampled_data, 
                             reverse = TRUE,
                             timewt = "n/G2")$concordance
    AUC_boot_1 <- roc((event)~pred_prob_1,ci=TRUE, data = resampled_data)
    
    
    harrell_C_boot_1 <- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                                    resampled_data, 
                                    reverse = TRUE)$concordance
    
    Uno_C_boot_1<- concordance(Surv(surv_time, event) ~ lp_bs_1, 
                               resampled_data, 
                               reverse = TRUE,
                               timewt = "n/G2")$concordance
    
    AUC_change <- AUC_boot_1$auc- AUC_boot$auc
    harrell_C_boot_change <- harrell_C_boot_1 -harrell_C_boot
    Uno_C_boot_change <- Uno_C_boot_1- Uno_C_boot
    #performance measures in the original data
    lp_orig<- predict(boot_md0, type = 'lp',newdata = impds)
    lp_orig_1 <- predict(boot_md1, type = 'lp',newdata = impds)
    pred_prob_vl <- predict(boot_md0, type = "risk", time = 365.25, newdata = impds)
    pred_prob_vl_1 <- predict(boot_md1, type = "risk", time = 365.25, newdata = impds)
    AUC_orig <- roc((event)~pred_prob_vl,ci=TRUE, data = impds)
    AUC_orig_1 <- roc((event)~pred_prob_vl_1,ci=TRUE, data = impds)
    harrell_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                                  impds, 
                                  reverse = TRUE)$concordance
    
    Uno_C_orig <- concordance(Surv(surv_time, event) ~ lp_orig, 
                              impds, 
                              reverse = TRUE,
                              timewt = "n/G2")$concordance
    harrell_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                    impds, 
                                    reverse = TRUE)$concordance
    
    Uno_C_orig_1 <- concordance(Surv(surv_time, event) ~ lp_orig_1, 
                                impds, 
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
    AUC_orig_change <- AUC_orig_1$auc - AUC_orig$auc
    harrell_C_orig_change <- harrell_C_orig_1- harrell_C_orig
    Uno_C_orig_change <- Uno_C_orig_1 - Uno_C_orig
    #optimism
    result[j,1] <- app_cstat_model_change - test_cstat_model_change
    result[j,2] <-  harrell_C_boot_change - harrell_C_orig_change
    result[j,3] <-   Uno_C_boot_change - Uno_C_orig_change
    result[j,4] <- AUC_orig_change
    result[j,5] <- harrell_C_orig_change
    result[j,6] <- Uno_C_orig_change
  }
  boot_md00 <- coxph(Surv(surv_time, event) ~  histology_coded + ps + 
                       NLR + firstline + age, data = impds)
  boot_md10 <- update(boot_md00, . ~ . + eq_5d)
  lp_app <- predict(boot_md00)
  lp_app_1 <- predict(boot_md10)
  pred_prob <- predict(boot_md00,type = "risk")
  pred_prob_1 <- predict(boot_md10,type = "risk")
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
  Uno_C0_1 <- concordance(Surv(surv_time, event) ~ lp_app_1, 
                          impds, 
                          reverse = TRUE,
                          timewt = "n/G2")$concordance
  
  app_uno_upci_1 <- Uno_C0_1 + 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  app_uno_loci_1 <- Uno_C0_1 - 1.96*sqrt(concordance(Surv(surv_time, event) ~ lp_app_1, 
                                                     impds, 
                                                     reverse = TRUE,
                                                     timewt = "n/G2")$var)
  
  Uno_C0_change <- Uno_C0_1 - Uno_C0
  harrell_C0_change <- concordance(boot_md10)$concordance -  concordance(boot_md00)$concordance
  auc0_change<- roc((surv_time<365.25)~pred_prob_1,ci=TRUE, data = impds)$auc - roc((surv_time<365.25)~pred_prob,ci=TRUE, data = impds)$auc
  
  #optimism corrected
  result_imp_c_1[i,1] <- auc0_change-mean(result[,1])
  result_imp_c_1[i,2:3] <- quantile(result[,4], c(0.025,0.975))
  result_imp_c_1[i,4] <- harrell_C0_change - mean(result[,2])
  result_imp_c_1[i,5:6] <- quantile(result[,5], c(0.025,0.975))
  result_imp_c_1[i,7] <- Uno_C0_change - mean(result[,3])
  result_imp_c_1[i,8:9] <- quantile(result[,6], c(0.025,0.975))
}
time2 <- Sys.time()
time2-time1
colnames(result_imp_c_1) <- c("AUC","AUC_lo","AUC_up","Harrell C", "HC_lo", "HC_up", "UnoC","UnoC_lo","UnoC_up")
summary(result_imp_c_1)
