rm(list = ls())
cat("\014")  

library(lme4)
library(car)

# Read data
an_data <- read.table("data_tables/Fig1_AN_data_all.txt", header = TRUE)

an_data$chinID <- as.factor(an_data$chinID)
an_data$hearing <- as.factor(an_data$hearing)
str(an_data)

## Model 1 
m_ttr <- lmer(TTR_dB  ~ cf_Hz_log * hearing + (1|chinID), data=an_data)
Anova(m_ttr, test.statistic='F')
