rm(list = ls())
cat("\014")  

restrict_cf <- function(data,minVal,maxVal){  
  data <- data[data$CF_kHz >= minVal,]
  data <- data[data$CF_kHz <= maxVal,]
  return(data)
}

restrict_sr <- function(data,minVal,maxVal){  
  data <- data[data$SpontRate >= minVal,]
  data <- data[data$SpontRate <= maxVal,]
  return(data)
}

library(lme4)
library(car)
dat_all <- read.table('../tables_for_stats/Fig8_table_fric_in_noise_DR.txt', header=TRUE)
dat_all$StimSNR <- as.ordered(dat_all$StimSNR)
dat_all$HearingStatus<- as.factor(dat_all$HearingStatus)
dat_all$CF_kHz_log <- log(dat_all$CF_kHz)
dat_all$SRgroup <- as.factor(dat_all$SpontRate > 18)

str(dat_all)



dat_hf <- restrict_cf(dat_all, 2.5, 8) 
dat_lsr <- restrict_sr(dat_all, 0, 18) 
dat_hsr <- restrict_sr(dat_all, 18, Inf) 


m_sn_all <- lm(DR_SN  ~ CF_kHz_log * HearingStatus + SRgroup * StimSNR * HearingStatus, data=dat_hf)
Anova(m_sn_all, test.statistic='F')

m_sn_lsr <- lm(DR_SN  ~ CF_kHz_log * HearingStatus + StimSNR * HearingStatus, data=dat_lsr)
Anova(m_sn_lsr, test.statistic='F')

m_sn_hsr <- lm(DR_SN  ~ CF_kHz_log * HearingStatus + StimSNR * HearingStatus, data=dat_hsr)
Anova(m_sn_hsr, test.statistic='F')



m_n_all <- lm(DR_N  ~ CF_kHz_log * HearingStatus + SRgroup * StimSNR * HearingStatus, data=dat_hf)
Anova(m_n_all, test.statistic='F')

m_n_lsr <- lm(DR_N  ~ CF_kHz_log * HearingStatus + StimSNR * HearingStatus, data=dat_lsr)
Anova(m_n_lsr, test.statistic='F')

m_sn_hsr <- lm(DR_N  ~ CF_kHz_log * HearingStatus + StimSNR * HearingStatus, data=dat_hsr)
Anova(m_sn_hsr, test.statistic='F')



