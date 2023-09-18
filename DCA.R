#Load packages
library(survival)
library(rms)
library(dplyr)
source("stdca.R")

#DCA in the internal validation data--------
#get the bootstrapping corrected predicted risk at 1 year
#1 hour to run
B <- 500
m <- 40
pred_md1 <- matrix(nrow = nrow(impds))
pred_md2 <- matrix(nrow = nrow(impds))
pred_md3 <- matrix(nrow = nrow(impds))
pred_md4 <- matrix(nrow = nrow(impds))

pred_iv_1 <- matrix(nrow = nrow(impds))
pred_iv_2 <- matrix(nrow = nrow(impds))
pred_iv_3 <- matrix(nrow = nrow(impds))
pred_iv_4 <- matrix(nrow = nrow(impds))
time1 <- Sys.time()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  
  set.seed(678)
  for (j in 1:B) {
    samp_index <- sample(1:nrow(impds), nrow(impds), rep=TRUE)
    resampled_data <- impds[samp_index,]
    #boot_model
    Non_eprom_model<- coxph(Surv(surv_time, event) ~ histology_coded + ps + 
                       NLR + firstline + age, data = resampled_data,  x = T)
    eq_5d_model <- update(Non_eprom_model,.~.+eq_5d)
    symp_burden_model <- update(Non_eprom_model,.~.+symburden)
    symp_count_model <- update(Non_eprom_model,.~.+symcount)
    
    pred.orig <- riskRegression::predictRisk(Non_eprom_model, 
                                                            newdata = impds,
                                                            times = 365.25)
    pred_eq <- riskRegression::predictRisk(eq_5d_model,
                                                  newdata = impds, 
                                                  times = 365.25)
    pred_symb <- riskRegression::predictRisk(symp_burden_model,
                                                    newdata = impds, 
                                                    times = 365.25)
    pred_symc <- riskRegression::predictRisk(symp_count_model,
                                                    newdata = impds, 
                                                    times = 365.25)
    
    pred_iv_1 <- cbind(pred_iv_1,pred.orig)
    pred_iv_2 <- cbind(pred_iv_2,pred_eq)
    pred_iv_3 <- cbind(pred_iv_3,pred_symb)
    pred_iv_4 <- cbind(pred_iv_4,pred_symc)
  }
  pred_md1  <- cbind(rowMeans(pred_iv_1[,-1]),pred_md1)
  pred_md2  <- cbind(rowMeans(pred_iv_2[,-1]),pred_md2)
  pred_md3  <- cbind(rowMeans(pred_iv_3[,-1]),pred_md3)
  pred_md4  <- cbind(rowMeans(pred_iv_4[,-1]),pred_md4)
}
time2 <- Sys.time()
time2-time1#1hour


##prepare the data to plot the DCA curve-----------

impds <- implong_pts_predict%>%filter(.imp == 0)

impds$pred1 <- rowMeans(pred_md1[,-(m+1)])
impds$pred1_eq <- rowMeans(pred_md2[,-(m+1)])
impds$pred1_symb <- rowMeans(pred_md3[,-(m+1)])
impds$pred1_symc <- rowMeans(pred_md4[,-(m+1)])
impds <- as.data.frame(impds)
dca_noneprom <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth <- smooth(dca_noneprom$net.benefit$pred1
                           [!is.na(dca_noneprom$net.benefit$pred1)],
                           twiceit = TRUE)
dca_smooth <- c(dca_smooth, 
                      rep(NA, sum(is.na(dca_noneprom$net.benefit$pred1))))



impds <- as.data.frame(impds)
dca_eq5d <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_eq", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_eq <- smooth(dca_eq5d$net.benefit$pred1_eq
                     [!is.na(dca_eq5d$net.benefit$pred1_eq)],
                     twiceit = TRUE)
dca_smooth_eq <- c(dca_smooth_eq, 
                rep(NA, sum(is.na(dca_eq5d$net.benefit$pred1_eq))))

dca_symb <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_symb", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_symb <- smooth(dca_symb$net.benefit$pred1_symb
                        [!is.na(dca_symb$net.benefit$pred1_symb)],
                        twiceit = TRUE)
dca_smooth_symb <- c(dca_smooth_symb, 
                   rep(NA, sum(is.na(dca_symb$net.benefit$pred1_symb))))

dca_symc <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_symc", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_symc <- smooth(dca_symc$net.benefit$pred1_symc
                          [!is.na(dca_symc$net.benefit$pred1_symc)],
                          twiceit = TRUE)
dca_smooth_symc <- c(dca_smooth_symc, 
                     rep(NA, sum(is.na(dca_symc$net.benefit$pred1_symc))))

#plot
png(filename = "DCA_vali.png", width = 3546, height = 2652, res = 400)
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_noneprom$net.benefit$threshold,
     dca_smooth,
     type = "l", 
     lwd = 3, 
     lty = 2,
     xlab = "Threshold probability", 
     ylab = "Net Benefit",
     xlim = c(0, 1),
     ylim = c(-0.10, 0.45), 
     bty = "n",
     cex.lab = 1.2, 
     cex.axis = 1,
     col = 4
)
abline(h = 0, 
       type = "l", 
       lwd = 3, 
       lty = 4,
       col = 8)
lines(dca_noneprom$net.benefit$threshold, 
      dca_noneprom$net.benefit$all, 
      type = "l",
      lwd = 3,
      lty = 5,
      col = 2)
lines(dca_eq5d$net.benefit$threshold, 
      dca_smooth_eq, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 7)
lines(dca_symb$net.benefit$threshold, 
      dca_smooth_symb, 
       type = "l", 
       lwd = 3, 
       lty = 5,
       col = 3)
lines(dca_symc$net.benefit$threshold, 
      dca_smooth_symc, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 6)
legend("topright",
       c(
         "Treat None",
         "Non-ePROM model",
         "EQ-5D model",
         "Symptom burden model",
         "Moderate-to-severe symptom model",
         "Treat All"
       ),
       lty = c(4, 5, 5, 5, 5, 4), 
       lwd = 3, 
       col = c(2, 4, 7, 3, 6, 8),
       bty = "n"
)
title("Decision curve analysis in the internal validation dataset", cex = 1.5)
dev.off()



#apparent -----
m <- 40
pred_app_1 <- matrix(nrow = nrow(impds),ncol = m)
pred_app_2 <- matrix(nrow = nrow(impds),ncol = m)
pred_app_3 <- matrix(nrow = nrow(impds),ncol = m)
pred_app_4 <- matrix(nrow = nrow(impds),ncol = m)
optimism_cali <- data.frame()
for(i in 1:m){
  impds <- implong_pts_predict%>%filter(.imp == i)
  md0 <-update(non_ePROM_md[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  md1 <-update(eprom_md1[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  md2 <-update(eprom_md2[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  md3 <-update(eprom_md3[["analyses"]][[i]], x=TRUE, y=TRUE, data = impds)
  impds$pred_app_1 <- riskRegression::predictRisk(md0, 
                                                newdata = impds,
                                                times = 365.25)
  impds$pred_app_2 <- riskRegression::predictRisk(md1, 
                                                newdata = impds,
                                                times = 365.25)
  impds$pred_app_3 <- riskRegression::predictRisk(md2, 
                                                newdata = impds,
                                                times = 365.25)
  impds$pred_app_4 <- riskRegression::predictRisk(md3, 
                                                newdata = impds,
                                                times = 365.25)
  
  pred_app_1[,i]  <- impds$pred_app_1
  pred_app_2[,i]  <- impds$pred_app_2
  pred_app_3[,i]  <- impds$pred_app_3
  pred_app_4[,i]  <- impds$pred_app_4
}



######prepare the data to plot---------
impds <- implong_pts_predict%>%filter(.imp == 0)

impds$pred1_app <- rowMeans(pred_app_1)
impds$pred1_eq_app <- rowMeans(pred_app_2)
impds$pred1_symb_app <- rowMeans(pred_app_3)
impds$pred1_symc_app <- rowMeans(pred_app_4)
impds <- as.data.frame(impds)
dca_noneprom_app <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_app", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_app <- smooth(dca_noneprom_app$net.benefit$pred1
                     [!is.na(dca_noneprom_app$net.benefit$pred1)],
                     twiceit = TRUE)
dca_smooth_app <- c(dca_smooth_app, 
                rep(NA, sum(is.na(dca_noneprom_app$net.benefit$pred1))))



impds <- as.data.frame(impds)
dca_eq5d_app <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_eq_app", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_eq_app <- smooth(dca_eq5d_app$net.benefit$pred1_eq
                        [!is.na(dca_eq5d_app$net.benefit$pred1_eq)],
                        twiceit = TRUE)
dca_smooth_eq_app <- c(dca_smooth_eq, 
                   rep(NA, sum(is.na(dca_eq5d_app$net.benefit$pred1_eq))))

dca_symb_app <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_symb_app", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_symb_app <- smooth(dca_symb_app$net.benefit$pred1_symb
                          [!is.na(dca_symb_app$net.benefit$pred1_symb)],
                          twiceit = TRUE)
dca_smooth_symb_app <- c(dca_smooth_symb, 
                     rep(NA, sum(is.na(dca_symb_app$net.benefit$pred1_symb))))

dca_symc_app <- stdca(
  data = impds, 
  outcome = "event", 
  ttoutcome = "surv_time",
  timepoint = 365.25, 
  predictors = "pred1_symc_app", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)
dca_smooth_symc_app <- smooth(dca_symc_app$net.benefit$pred1_symc
                          [!is.na(dca_symc_app$net.benefit$pred1_symc)],
                          twiceit = TRUE)
dca_smooth_symc_app <- c(dca_smooth_symc, 
                     rep(NA, sum(is.na(dca_symc_app$net.benefit$pred1_symc))))


png(filename = "DCA_dev.png", width = 3546, height = 2652, res = 400)
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_noneprom_app$net.benefit$threshold,
     dca_smooth_app,
     type = "l", 
     lwd = 3, 
     lty = 2,
     xlab = "Threshold probability", 
     ylab = "Net Benefit",
     xlim = c(0, 1),
     ylim = c(-0.10, 0.45), 
     bty = "n",
     cex.lab = 1.2, 
     cex.axis = 1,
     col = 4
)
abline(h = 0, 
       type = "l", 
       lwd = 3, 
       lty = 4,
       col = 8)
lines(dca_noneprom_app$net.benefit$threshold, 
      dca_noneprom_app$net.benefit$all, 
      type = "l",
      lwd = 3,
      lty = 5,
      col = 2)
lines(dca_eq5d_app$net.benefit$threshold, 
      dca_smooth_eq_app, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 7)
lines(dca_symb_app$net.benefit$threshold, 
      dca_smooth_symb_app, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 3)
lines(dca_symc_app$net.benefit$threshold, 
      dca_smooth_symc_app, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 6)
legend("topright",
       c(
         "Treat None",
         "Non-ePROM model",
         "EQ-5D model",
         "Symptom burden model",
         "Moderate-to-severe symptom model",
         "Treat All"
       ),
       lty = c(4, 5, 5, 5, 5, 4), 
       lwd = 3, 
       col = c(2, 4, 7, 3, 6, 8),
       bty = "n"
)
title("Decision curve analysis in the derivation dataset", cex = 1.5)
dev.off()



cbind(dca_noneprom$net.benefit,dca_eq5d$net.benefit)
