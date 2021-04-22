rm(list = ls())
cat("\014")  

library(lme4)
library(car)

# Read data
dp_data <- read.table("data_tables/Fig1_AN_data_all.txt", header = TRUE)

dp_data$chinID <- as.factor(dp_data$chinID)
dp_data$hearing <- as.factor(dp_data$hearing)
str(dp_data)

## Model 1 
m_thresh <- lm(an_thresh  ~ cf_Hz_log * hearing, data=dp_data)
Anova(m_thresh, test.statistic='F')

m_q10 <- lm(Q10local ~ cf_Hz_log * hearing, data=dp_data)
Anova(m_q10, test.statistic='F')
