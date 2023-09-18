#imputation
library(mice)
library(dplyr)

predict_ds <- survSplit(Surv(surv_time, event) ~ ., data = pts_predict, cut = 365.25,
episode="timegroup")

pts_predict_1 <- subset(predict_ds, timegroup == 1) 
pts_predict_1 <- cbind(pts_predict_1,pts_predict[,c("surv_time","event")])
colnames(pts_predict_1)[17] <- "surv_time_obs"
colnames(pts_predict_1)[18] <- "event_obs"
pts_predict
str(pts_all)
pts_all$sex<- as.numeric(pts_all$sex)-1
pts_all$SmokingHistory <- as.factor(pts_all$SmokingHistory)
pts_all$histology_coded <- as.factor(pts_all$histology_coded)
table(pts_predict$histology_coded)
pts_predict <- pts_all %>% select(-c(ChristieNo))
str(pts_predict)

ini <- mice(pts_predict_1,m=1, maxit = 464, seed = 1234, remove.constant = FALSE)
pred_matrix <- ini$predictorMatrix
pred_matrix[c("event","surv_time"),] <- 0

time1 <- Sys.time()
imp_pts_predict <- mice(pts_predict_1,m=40, maxit = 464, seed = 1234, predictorMatrix = pred_matrix,remove.constant = FALSE)
time2 <- Sys.time()
timediff <- time2 - time1
timediff
implong_pts_predict <- complete(imp_pts_predict,"long",include=TRUE)

#diagnostic plot

densityplot(imp_pts_predict,~  PDL1Score + NLR +eq_vas +eq_5d +symburden + symcount, thicker = 5)


