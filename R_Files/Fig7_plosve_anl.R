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
data_an <- read.table("data_tables/Fig7_AN_plosive_all.txt", header = TRUE)
data_an$ChinID <- as.factor(data_an$ChinID)
data_an$HearingStatus<- as.factor(data_an$HearingStatus)
data_an$SRgroup <- as.factor(data_an$SpontRate > 18)
data_an$CF_kHz_log <- log(data_an$CF_kHz)
str(data_an)

data_g <- restrict(data_an, 0.5, 3)
data_d <- restrict(data_an, 0.5, 8)

print("Start: ------------- The effect of hearing loss") 
### Without SpontRate
## (/g/) 
# Onset rate 
m_OR_g <- lm(OnRate_g ~ HearingStatus*CF_kHz_log, data=data_g)
Anova(m_OR_g, test.statistic='F')


# Driven rate 
m_DR_g <- lm(SusRate_g ~ HearingStatus*CF_kHz_log, data=data_g)
Anova(m_DR_g, test.statistic='F')


## (/d/) 
# Onset rate 
m_OR_d <- lm(OnRate_d ~ HearingStatus*CF_kHz_log, data=data_d)
Anova(m_OR_d, test.statistic='F')


# Driven rate 
m_DR_d <- lm(SusRate_d ~ HearingStatus*CF_kHz_log, data=data_d)
Anova(m_DR_d, test.statistic='F')

print("End: ------------- The effect of hearing loss") 



# ### With SpontRate
# print("Start : ------------- Factors contributing to coding deficit") 
# ## (/g/) 
# # Onset rate 
# m_OR_g <- lmer(OnRate_g ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_g)
# Anova(m_OR_g, test.statistic='F')
# 
# # Driven rate 
# m_DR_g <- lmer(SusRate_g ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_g)
# Anova(m_DR_g, test.statistic='F')
# 
# ## (/d/) 
# # Onset rate 
# m_OR_d <- lmer(OnRate_d ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_d)
# Anova(m_OR_d, test.statistic='F')
# 
# # Driven rate 
# m_DR_d <- lmer(SusRate_d ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_d)
# Anova(m_DR_d, test.statistic='F')


# ### With SpontRate (full CF range )
# ## (/g/) 
# # Onset rate 
# m_OR_g <- lmer(OnRate_g ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_an)
# Anova(m_OR_g, test.statistic='F')
# 
# # Driven rate 
# m_DR_g <- lmer(SusRate_g ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_an)
# Anova(m_DR_g, test.statistic='F')
# 
# ## (/d/) 
# # Onset rate 
# m_OR_d <- lmer(OnRate_d ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_an)
# Anova(m_OR_d, test.statistic='F')
# 
# # Driven rate 
# m_DR_d <- lmer(SusRate_d ~ HearingStatus*CF_kHz_log + SpontRate + (1|ChinID), data=data_an)
# Anova(m_DR_d, test.statistic='F')


# ### Without SpontRate (full CF range )
# ## (/g/) 
# # Onset rate 
# m_OR_g <- lmer(OnRate_g ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_an)
# Anova(m_OR_g, test.statistic='F')
# 
# # Driven rate 
# m_DR_g <- lmer(SusRate_g ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_an)
# Anova(m_DR_g, test.statistic='F')
# 
# ## (/d/) 
# # Onset rate 
# m_OR_d <- lmer(OnRate_d ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_an)
# Anova(m_OR_d, test.statistic='F')
# 
# # Driven rate 
# m_DR_d <- lmer(SusRate_d ~ HearingStatus*CF_kHz_log + (1|ChinID), data=data_an)
# Anova(m_DR_d, test.statistic='F')


## Driven Rate re. Spont 
data_an$OnRate_g_rel <- data_an$OnRate_g - data_an$SpontRate
data_an$SusRate_g_rel <- data_an$SusRate_g - data_an$SpontRate
data_an$OnRate_d_rel <- data_an$OnRate_d - data_an$SpontRate
data_an$SusRate_d_rel <- data_an$SusRate_d - data_an$SpontRate

data_g$OnRate_g_rel <- data_g$OnRate_g - data_g$SpontRate
data_g$SusRate_g_rel <- data_g$SusRate_g - data_g$SpontRate
data_g$OnRate_d_rel <- data_g$OnRate_d - data_g$SpontRate
data_g$SusRate_d_rel <- data_g$SusRate_d - data_g$SpontRate

data_d$OnRate_g_rel <- data_d$OnRate_g - data_d$SpontRate
data_d$SusRate_g_rel <- data_d$SusRate_g - data_d$SpontRate
data_d$OnRate_d_rel <- data_d$OnRate_d - data_d$SpontRate
data_d$SusRate_d_rel <- data_d$SusRate_d - data_d$SpontRate


## (/g/) 
# Onset rate 
m_OR_g <- lm(OnRate_g_rel ~ HearingStatus*CF_kHz_log, data=data_g)
Anova(m_OR_g, test.statistic='F')

# Driven rate 
m_DR_g <- lm(SusRate_g_rel ~ HearingStatus*CF_kHz_log, data=data_g)
Anova(m_DR_g, test.statistic='F')

## (/d/) 
# Onset rate 
m_OR_d <- lm(OnRate_d_rel ~ HearingStatus*CF_kHz_log, data=data_d)
Anova(m_OR_d, test.statistic='F')

# Driven rate
m_DR_d <- lm(SusRate_d_rel ~ HearingStatus*CF_kHz_log, data=data_d)
Anova(m_DR_d, test.statistic='F')

####################### FFR 
# Read data
data_ffr <- read.table("data_tables/Fig7_FFR_plosive_all.txt", header = TRUE)
data_ffr$HearingStatus <- as.factor(data_ffr$HearingStatus)
data_ffr$ChinID <- as.factor(data_ffr$ChinID)
str(data_ffr)

m_FFR_d <- lmer(OnRate_d ~ HearingStatus + (1|ChinID), data=data_ffr)
Anova(m_FFR_d, test.statistic='F')

m_FFR_g <- lmer(OnRate_g ~ HearingStatus + (1|ChinID), data=data_ffr)
Anova(m_FFR_g, test.statistic='F')
