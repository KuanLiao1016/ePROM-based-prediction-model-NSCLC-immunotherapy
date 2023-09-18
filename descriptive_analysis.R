#Descriptive analysis
library(dplyr)
library(ggplot2)
library(mice)
library(survminer)
pts_all%>%
  summarise(
    mean(age,na.rm = T),
    sd(age,na.rm = T),
    quantile(PDL1Score,c(0.25,0.5,0.75),na.rm = T),
    mean(NLR,na.rm = T),
    sd(NLR,na.rm = T),
    VAS = quantile(eq_vas,c(0.25,0.5,0.75),na.rm = T),
    EQ5D = quantile(eq_5d,c(0.25,0.5,0.75),na.rm = T),
    symptom = quantile(symburden, c(0.25,0.5,0.75),na.rm = T),
    symptomcount = quantile(symcount, c(0.25,0.5,0.75),na.rm = T),
  )
md.pattern(pts_all)
#sex
pts_all <- pts_all%>%mutate(
  sex = as.factor(case_when(sex == 'Female'~ 0,
                  sex == 'Male' ~ 1)))

pts_all%>%group_by(sex)%>%
  summarise(n=n())%>%
  mutate(freq =  round(n/sum(n)*100,1))

#ECOG PS
pts_all%>%group_by(ps)%>%
  summarise(n=n())%>%
  mutate(freq =  round(n/sum(n)*100,1))


#PDL1
pts_all <- pts_all%>%mutate(
  PDL1group = case_when(PDL1Score<1 ~ 1,
                        PDL1Score>=1&PDL1Score<50 ~ 2,
                        PDL1Score>=50 ~ 3))

pts_all%>%group_by(PDL1group)%>%
  summarise(n=n())%>%
  mutate(freq = round(n/439*100,1))

#NLR
pts_all%>%group_by(is.na(NLR))%>%
  summarise(n=n())%>%
  mutate(freq =  round(n/sum(n)*100,2))


#histology
pts_all%>%group_by(histology_coded)%>%
  summarise(n=n())%>%
  mutate(freq = round(n/sum(n)*100,1))
#smoking history
pts_all%>%group_by(SmokingHistory)%>%
  summarise(n=n())%>%
  mutate(freq =   round(n/463*100,1))
#line of treatment
pts_all%>%group_by(firstline)%>%
  summarise(n=n())%>%
  mutate(freq =   round(n/sum(n)*100,1))



#moderate to severe symptoms
pts_all <- pts_all%>%mutate(
  symcountgroup = case_when(symcount==0 ~ 0,
                            symcount==1 ~ 1,
                        symcount>=2 ~ 2))

pts_all%>%group_by(symcountgroup)%>%
  summarise(n=n(),)%>%
  mutate(freq = round(n/186*100,1))




#survival outcome KM

surv_outcome <- survfit(Surv(pts_all$surv_time/30.4167,pts_all$event) ~ 1, data = pts_all)
surv_plot <- ggsurvplot(surv_outcome,
                          pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, # Add risk table
                          linetype = "strata", # Change line type by groups
                          surv.median.line = "hv", # Specify median survival
                          ggtheme = theme_bw(), # Change ggplot2 theme
                          palette = c("#454545"),
                          ylab = 'Overall survival',
                          xlab = 'Months',
                          legend='none')

time_of_interest <- 12

df <- data.frame(
  times = time_of_interest,
  probs = summary(surv_outcome, time_of_interest)$surv
)

surv_plot$plot +
  geom_segment(data = df,
               aes(x = times, y = 0, xend = times, yend = probs), linetype = "dashed", colour = "red")+
  geom_segment(data = df,
               aes(x = 0, y = probs, xend = times, yend = probs), linetype = "dashed", colour = "red")
ggsave("KM.png",dpi = 500)
summary(surv_outcome, times = 12)

ggsave("km_plot.png",dpi = 500)
