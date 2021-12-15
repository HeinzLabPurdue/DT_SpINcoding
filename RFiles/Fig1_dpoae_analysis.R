rm(list = ls())
cat("\014")  

library(lme4)
library(car)

# Read data
dp_data <- read.table("../tables_for_stats/Fig1_dpoae_data_all.txt", header = TRUE)

dp_data$chinID <- as.factor(dp_data$chinID)
dp_data$hearing <- as.factor(dp_data$hearing)
dp_data$freq_kHz <- as.factor(dp_data$freq_kHz)
str(dp_data)

## Model 1 
m_DPOAE <- lmer(dpoae_amp  ~ freq_kHz * hearing + (1|chinID), data=dp_data)
Anova(m_DPOAE, test.statistic='F')
# summary(m_DPOAE)
# anova(m_DPOAE)