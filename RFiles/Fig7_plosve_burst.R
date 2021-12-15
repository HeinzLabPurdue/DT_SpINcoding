rm(list = ls())
cat("\014")  

restrict <- function(f,minVal,maxVal){  
  f <- f[f$CF_kHz >= minVal,]
  f <- f[f$CF_kHz <= maxVal,]
  return(f)
}

library(lme4)
library(car)

# Read data
data_an <- read.table("../tables_for_stats/Fig7_AN_plosive_burst.txt", header = TRUE)
data_an$ChinID <- as.factor(data_an$ChinID)
data_an$HearingStatus<- as.factor(data_an$HearingStatus)
data_an$SRgroup <- as.factor(data_an$SpontRate > 18)
data_an$CF_kHz_log <- log(data_an$CF_kHz)
str(data_an)

data_d <- restrict(data_an, 0.5, 8)
data_g <- restrict(data_an, 0.5, 3)

print("Start: ------------- The effect of hearing loss") 
### Without SpontRate
# m_DR_g <- lmer(SusRate_g ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_g)
# Anova(m_DR_g, test.statistic='F')

## (/d/) 

# Driven rate 
m_DR_d <- lm(BurstRates_d  ~ HearingStatus*CF_kHz_log, data=data_d)
Anova(m_DR_d, test.statistic='F')

## (/g/) 

# Driven rate 
m_DR_g <- lm(BurstRates_g ~ HearingStatus*CF_kHz_log, data=data_g)
Anova(m_DR_g, test.statistic='F')


# m_DR_d <- lmer(SusRate_d ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_d)
# Anova(m_DR_d, test.statistic='F')
print("End: ------------- The effect of hearing loss") 



####################### FFR 
# Read data
data_ffr <- read.table("../tables_for_stats/Fig7_FFR_plosive_all.txt", header = TRUE)
data_ffr$HearingStatus <- as.factor(data_ffr$HearingStatus)
data_ffr$ChinID <- as.factor(data_ffr$ChinID)
str(data_ffr)

m_FFR_d <- lmer(OnRate_d ~ HearingStatus + (1|ChinID), data=data_ffr)
Anova(m_FFR_d, test.statistic='F')

m_FFR_g <- lmer(OnRate_g ~ HearingStatus + (1|ChinID), data=data_ffr)
Anova(m_FFR_g, test.statistic='F')
