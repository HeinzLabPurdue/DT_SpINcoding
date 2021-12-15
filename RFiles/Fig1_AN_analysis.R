rm(list = ls())
cat("\014")  

library(lme4)
library(car)

# Read data
dp_data <- read.table("../tables_for_stats/Fig1_AN_data_all.txt", header = TRUE)

dp_data$chinID <- as.factor(dp_data$chinID)
dp_data$hearing <- as.factor(dp_data$hearing)
str(dp_data)

## Model 1 
m_thresh <- lm(an_thresh  ~ cf_Hz_log * hearing, data=dp_data)
Anova(m_thresh, test.statistic='F')
# summary(m_thresh)
# anova(m_thresh)

m_q10 <- lm(Q10local ~ cf_Hz_log * hearing, data=dp_data)
Anova(m_q10, test.statistic='F')

# 
# ## Model 1 
# m_thresh <- lmer(an_thresh  ~ cf_Hz_log * hearing + (1|chinID), data=dp_data)
# Anova(m_thresh, test.statistic='F')
# # summary(m_thresh)
# # anova(m_thresh)
# 
# m_q10 <- lmer(Q10local ~ cf_Hz_log * hearing + (1|chinID), data=dp_data)
# Anova(m_q10, test.statistic='F')
