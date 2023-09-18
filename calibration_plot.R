#library packages
library(ggplot2)
library(rms)
library(survival)
library(dplyr)
library(mice)
library(riskRegression)
library(cowplot)
#---- Non-ePROM model ----

#optimism corrected predicted risk
B <- 500
m <- 40
pred <- matrix(nrow = nrow(impds))
inter_cali <- matrix(nrow = nrow(impds))
time1 <- Sys.time()
for(i in 1:m){
impds <- implong_pts_predict%>%filter(.imp == i)

set.seed(678)
for (j in 1:B) {
  samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
  resampled_data <- impds[samp_index,]
  #boot_model
  boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                     NLR + firstline + age, data = resampled_data,  x = T)
  

  resampled_data$pred.orig <- riskRegression::predictRisk(boot_md0, 
                                                     newdata = impds,
                                                     times = 365.25)
  
  inter_cali<- cbind(inter_cali,resampled_data$pred.orig)
}
pred  <- cbind(rowMeans(inter_cali[,-1]),pred)
}
time2 <- Sys.time()
time2-time1

pred[,-(m+1)]

impds$pred.cll <- log(-log(1-rowMeans(pred[,-(m+1)])))
# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred[,-(m+1)])
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

#calibration curve- internal validation non-ePROM model
A_iv <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag = "A",
       title = "Non-ePROM model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))

A_iv_x <- cowplot::axis_canvas(A_iv, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

A_iv_combined <- insert_xaxis_grob(
  A_iv,
  A_iv_x,
  position = "bottom"
)
A_iv_combined <- ggdraw(A_iv_combined)
A_iv_combined


#apparent calibration plot - non-ePROM model
m <- 40
pred_app <- matrix(nrow = nrow(impds),ncol = m)
optimism_cali <- data.frame()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  cali_md1 <-update(non_ePROM_md[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  impds$pred_app <- riskRegression::predictRisk(cali_md1, 
                                            newdata = impds,
                                            times = 365.25)
 
  pred_app[,i]  <- impds$pred_app
}
impds$pred_app.cll <- log(-log(1-rowMeans(pred_app)))

# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred_app.cll,3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred_app)
  )

A_app <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag = "A",
       title = "Non-ePROM model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))
A_app_x <- cowplot::axis_canvas(A_app, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

A_app_combined <- insert_xaxis_grob(
  A_app,
  A_app_x,
  position = "bottom"
)
A_app_combined <- ggdraw(A_app_combined)
A_app_combined

#---- EQ-5D model----
#optimism corrected predicted risk
B <- 500
m <- 40
pred <- matrix(nrow = nrow(impds))
inter_cali <- matrix(nrow = nrow(impds))
time1 <- Sys.time()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  
  set.seed(678)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + eq_5d, data = resampled_data,  x = T)
    
    
    resampled_data$pred.orig <- riskRegression::predictRisk(boot_md0, 
                                                            newdata = impds,
                                                            times = 365.25)
    
    inter_cali<- cbind(inter_cali,resampled_data$pred.orig)
  }
  pred  <- cbind(rowMeans(inter_cali[,-1]),pred)
}
time2 <- Sys.time()
time2-time1


impds$pred.cll <- log(-log(1-rowMeans(pred[,-(m+1)])))
# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred[,-(m+1)])
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

#calibration curve- internal validation non-ePROM model
B_iv <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
    xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag= "B",
       title = "EQ-5D model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))


B_iv_x <- cowplot::axis_canvas(B_iv, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

B_iv_combined <- insert_xaxis_grob(
  B_iv,
  B_iv_x,
  position = "bottom"
)
B_iv_combined <- ggdraw(B_iv_combined)
B_iv_combined


#apparent calibration plot - EQ-5D model
m <- 40
pred_app <- matrix(nrow = nrow(impds),ncol = m)
optimism_cali <- data.frame()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  cali_md1 <-update(eprom_md1[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  impds$pred_app <- riskRegression::predictRisk(cali_md1, 
                                                newdata = impds,
                                                times = 365.25)
  
  pred_app[,i]  <- impds$pred_app
}
impds$pred_app.cll <- log(-log(1-rowMeans(pred_app)))

# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred_app.cll,3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred_app)
)

B_app <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag= "B",
       title = "EQ-5D model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))

B_app_x <- cowplot::axis_canvas(B_app, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

B_app_combined <- insert_xaxis_grob(
  B_app,
  B_app_x,
  position = "bottom"
)
B_app_combined <- ggdraw(B_app_combined)
B_app_combined


#---- Symptom burden model----
#optimism corrected predicted risk
B <- 500
m <- 40
pred <- matrix(nrow = nrow(impds))
inter_cali <- matrix(nrow = nrow(impds))
time1 <- Sys.time()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  
  set.seed(678)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + symburden, data = resampled_data,  x = T)
    
    
    resampled_data$pred.orig <- riskRegression::predictRisk(boot_md0, 
                                                            newdata = impds,
                                                            times = 365.25)
    
    inter_cali<- cbind(inter_cali,resampled_data$pred.orig)
  }
  pred  <- cbind(rowMeans(inter_cali[,-1]),pred)
}
time2 <- Sys.time()
time2-time1

impds$pred.cll <- log(-log(1-rowMeans(pred[,-(m+1)])))
# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred[,-(m+1)])
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

#calibration curve- internal validation Symptom burden model
C_iv <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       title = "Symptom burden model",
       tag = "C")+
  theme_minimal() +
  theme(legend.position = "right",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))

C_iv_x <- cowplot::axis_canvas(C_iv, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

C_iv_combined <- insert_xaxis_grob(
  C_iv,
  C_iv_x,
  position = "bottom"
)
C_iv_combined <- ggdraw(C_iv_combined)
C_iv_combined



#apparent calibration plot - Symptom burden model
m <- 40
pred_app <- matrix(nrow = nrow(impds),ncol = m)
optimism_cali <- data.frame()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  cali_md1 <-update(eprom_md2[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  impds$pred_app <- riskRegression::predictRisk(cali_md1, 
                                                newdata = impds,
                                                times = 365.25)
  
  pred_app[,i]  <- impds$pred_app
}
impds$pred_app.cll <- log(-log(1-rowMeans(pred_app)))

# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred_app.cll,3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred_app)
)

C_app <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       title = "Symptom burden model",
       tag = "C")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))


C_app_x <- cowplot::axis_canvas(C_app, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

C_app_combined <- insert_xaxis_grob(
  C_app,
  C_app_x,
  position = "bottom"
)
C_app_combined <- ggdraw(C_app_combined)
C_app_combined



#---- Moderate-to-severe symptom model----
#optimism corrected predicted risk
B <- 500
m <- 40
pred <- matrix(nrow = nrow(impds))
inter_cali <- matrix(nrow = nrow(impds))
time1 <- Sys.time()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  
  set.seed(678)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    resampled_data <- impds[samp_index,]
    #boot_model
    boot_md0<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age + symcount, data = resampled_data,  x = T)
    
    
    resampled_data$pred.orig <- riskRegression::predictRisk(boot_md0, 
                                                            newdata = impds,
                                                            times = 365.25)
    
    inter_cali<- cbind(inter_cali,resampled_data$pred.orig)
  }
  pred  <- cbind(rowMeans(inter_cali[,-1]),pred)
}
time2 <- Sys.time()
time2-time1

impds$pred.cll <- log(-log(1-rowMeans(pred[,-(m+1)])))
# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred[,-(m+1)])
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

#calibration curve- internal validation Symptom burden model
D_iv <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag= "D",
       title = "Moderate-to-severe symptom model") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))


D_iv_x <- cowplot::axis_canvas(D_iv, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

D_iv_combined <- insert_xaxis_grob(
  D_iv,
  D_iv_x,
  position = "bottom"
)
D_iv_combined <- ggdraw(D_iv_combined)
D_iv_combined

#apparent calibration plot - Symptom burden model
m <- 40
pred_app <- matrix(nrow = nrow(impds),ncol = m)
optimism_cali <- data.frame()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  cali_md1 <-update(eprom_md3[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  impds$pred_app <- riskRegression::predictRisk(cali_md1, 
                                                newdata = impds,
                                                times = 365.25)
  
  pred_app[,i]  <- impds$pred_app
}
impds$pred_app.cll <- log(-log(1-rowMeans(pred_app)))

# Estimate actual risk
vcal <- rms::cph(Surv(surv_time, event) ~ rcs(pred_app.cll,3),
                 x = T,
                 y = T,
                 surv = T,
                 data = impds
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal, 
                           times = 365.25, 
                           newdata = impds)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 365.25, 
                             newdata = impds)$lower,
  "pred" = rowMeans(pred_app)
)

D_app <- ggplot(dat_cal, aes(x = pred, y = obs)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey80") +
  geom_line(size = 0.5, linetype = "solid", color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "#d11141", size = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Predicted probability of 1-year mortality",
       y = "Observed probability of 1-year mortality",
       tag= "D",
       title = "Moderate-to-severe symptom model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8.5))

D_app_x <- cowplot::axis_canvas(D_app, axis = "x") +
  geom_histogram(
    data = dat_cal,
    aes(x = pred),
    colour = "grey50",
    bins = 35) 

D_app_combined <- insert_xaxis_grob(
  D_app,
  D_app_x,
  position = "bottom"
)
D_app_combined <- ggdraw(D_app_combined)
D_app_combined

#########Combine the plots#######
ggpubr::ggarrange(A_app_combined,B_app_combined,C_app_combined,D_app_combined,ncol = 2, nrow = 2)
ggsave("cali_plot1.png",dpi = 500)
ggpubr::ggarrange(A_iv_combined,B_iv_combined,C_iv_combined,D_iv_combined,ncol = 2, nrow = 2)
ggsave("cali_plot2.png",dpi = 500)
