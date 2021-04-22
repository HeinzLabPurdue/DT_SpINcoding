rm(list = ls())
cat("\014")  

library(lme4)
library(car)

# Read data
VPdat <- read.table("VPdistTab_50ms_1000ps.txt", header = TRUE)
VPdat$UnitID <- as.factor(VPdat$UnitID)
VPdat$HearingStatus <- as.factor(VPdat$HearingStatus)
VPdat$WindowDur_ms <- as.factor(VPdat$WindowDur_ms)
VPdat$VPcost <- as.factor(VPdat$VPcost)
str(VPdat)

## Model 1 
m_vp <- lm(VPdistance  ~ CF_Hz_log * HearingStatus * Rate, data=VPdat)
Anova(m_vp, test.statistic='F')
summary(m_vp)
plot(m_vp)