rm(list = ls())
cat("\014")  

get_deq <- function(mdl){  
  var_out <- data.frame(d_eq_ttr=0, d_eq_q10=0, d_eq_thr=0)
  
  depVars <- rownames(mdl)
  ttr_ind <- match("TTR_dB", depVars)
  Q10_ind <- match("Q10local", depVars)
  thr_ind <- match("Threshold_dBSPL", depVars)
  
  pVal_ttr= mdl$`Pr(>F)`[ttr_ind];
  df_ttr= mdl$Df.res[ttr_ind];
  var_out$d_eq_ttr = 2*qt(1 - pVal_ttr, df_ttr) / sqrt(df_ttr)
  
  pVal_q10= mdl$`Pr(>F)`[Q10_ind];
  df_q10= mdl$Df.res[Q10_ind];
  var_out$d_eq_q10 = 2*qt(1 - pVal_q10, df_q10) / sqrt(df_q10)
  
  pVal_thr= mdl$`Pr(>F)`[thr_ind];
  df_thr= mdl$Df.res[thr_ind];
  var_out$d_eq_thr = 2*qt(1 - pVal_thr, df_thr) / sqrt(df_thr)
  
  return(var_out)
}

library(lme4)
library(car)
dat <- read.table('../tables_for_stats/Fig8_table_fric_in_noise_DR.txt', header=TRUE)
str(dat)
dat$ChinID <- as.factor(dat$ChinID) # Value doesn't matter 
dat$UnitID <- as.factor(dat$UnitID) # Value doesn't matter 
dat$RespCorrTXF <- log(atanh(abs(dat$RespCorr))) 
dat$RespCorr_RawTXF <- log(1+atanh(.99*(dat$RespCorrRaw-min(dat$RespCorrRaw))/(max(dat$RespCorrRaw)-min(dat$RespCorrRaw)) ))
# dat$RespCorr_RawTXF <- atanh(dat$RespCorr)
# dat$RespCorr_RawTXF <- log(atanh(abs(dat$RespCorrRaw)))
dat$HearingStatus<- as.factor(dat$HearingStatus)
dat$SRgroup <- as.factor(dat$SpontRate > 18)
dat$StimSNR <- as.ordered(dat$StimSNR) # Value matters 
dat$CF_kHz_log <- log(dat$CF_kHz)
dat$SignifCoding <- as.factor(dat$RespCorrRaw > .001)

####### Start with raw and check floor at the end 

# Model #1
m4 <- lmer(RespCorrRaw  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus + (1|ChinID), data=dat)
Anova(m4, test.statistic='F')
plot(m4)

## Final
print("Final for Correlation")
m4 <- lm(RespCorr  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus, data=dat)
Anova(m4, test.statistic='F')


## Final
print("Final for DR")
mDR <- lm(DR_SNrelN  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus, data=dat)
Anova(mDR, test.statistic='F')


plot(m4)

# m4 <- lmer(RespCorr  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus + (1|ChinID), data=dat)
# Anova(m4, test.statistic='F')
# plot(m4)

# m4 <- lmer(RespCorr  ~ CF_kHz_log + StimSNR + SpontRate*HearingStatus + (1|ChinID), data=dat)
# Anova(m4, test.statistic='F')
# plot(m4)


# Model #2
m1_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + SRgroup + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m1_fricIN, test.statistic='F')
plot(m1_fricIN)

# Remove Threshold_dBSPL because lowest p-Value 
m2_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + SRgroup + Q10local + TTR_dB + (1|ChinID), data=dat)
stat_m2_fricIN <-Anova(m2_fricIN, test.statistic='F')
print(stat_m2_fricIN)
plot(m2_fricIN)

# Remove SRgroup because lowest p-Value 
m2_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Q10local + TTR_dB + (1|ChinID), data=dat)
stat_m2_fricIN <-Anova(m2_fricIN, test.statistic='F')
print(stat_m2_fricIN)

plot(m2_fricIN)

depVars <- rownames(stat_m2_fricIN)
ttr_ind <- match("TTR_dB", depVars)
Q10_ind <- match("Q10local", depVars)

pVal_ttr= stat_m2_fricIN$`Pr(>F)`[ttr_ind];
df_ttr= stat_m2_fricIN$Df.res[ttr_ind];
d_eq_ttr = 2*qt(1 - pVal_ttr, df_ttr) / sqrt(df_ttr)
r_eq_ttr = sqrt(qt(1 - pVal_ttr, df_ttr)^2 / (qt(1 - pVal_ttr, df_ttr)^2 + df_ttr^2))

pVal_q10= stat_m2_fricIN$`Pr(>F)`[Q10_ind];
df_q10= stat_m2_fricIN$Df.res[Q10_ind];
d_eq_q10 = 2*qt(1 - pVal_q10, df_q10) / sqrt(df_q10)
r_eq_q10 = sqrt(qt(1 - pVal_q10, df_q10)^2 / (qt(1 - pVal_q10, df_q10)^2 + df_q10^2))

sprintf("Effect size for [TTR= %.2f] [Q10= %.2f]", d_eq_ttr, d_eq_q10) 

# # Remove Q10 because lowest p-Value 
# m2_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + SRgroup + TTR_dB + (1|ChinID), data=dat)
# Anova(m2_fricIN, test.statistic='F')
# plot(m2_fricIN)

m2_FloorCheck <- lmer(RespCorr ~ CF_kHz_log + StimSNR + SRgroup + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m2_FloorCheck, test.statistic='F')
plot(m2_FloorCheck)
print("~similar effects for floored and unfloored data -----------------")

m2_FloorCheck <- lmer(RespCorr ~ CF_kHz_log + StimSNR + SpontRate + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m2_FloorCheck, test.statistic='F')
plot(m2_FloorCheck)
print("~similar effects for floored and unfloored data -----------------")


## deq check 
m1_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), data=dat)
anOut_f1_hf <- Anova(m1_fricIN, test.statistic='F')
print(anOut_f1_hf)
var_out <- get_deq(anOut_f1_hf)
print(var_out)

m1_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), data=dat)
anOut_f1_hf <- Anova(m1_fricIN, test.statistic='F')
print(anOut_f1_hf)
var_out <- get_deq(anOut_f1_hf)
print(var_out)

## Pruning model 
m1_reml <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), REML = FALSE, data=dat)
m2_reml <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), REML = FALSE, data=dat)
anova(m1_reml, m2_reml)


print("Ignore below this -----------------")
####### Start with floor and check raw at the end 
m4 <- lmer(RespCorr  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus + (1|ChinID), data=dat)
Anova(m4, test.statistic='F')
plot(m4)


m5 <- lmer(RespCorr ~ CF_kHz_log + StimSNR + SRgroup + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m5, test.statistic='F')
plot(m5)

m5 <- lmer(RespCorr ~ CF_kHz_log + StimSNR + SRgroup + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m5, test.statistic='F')
plot(m5)

m5 <- lmer(RespCorr ~ CF_kHz_log + StimSNR + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m5, test.statistic='F')
plot(m5)

m5 <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m5, test.statistic='F')
plot(m5)


##################
m5 <- lmer(RespCorrTXF ~ CF_kHz_log + StimSNR  + SRgroup + Q10local + TTR_dB + (1|ChinID), data=dat)
Anova(m5, test.statistic='F')
plot(m5)

# 
# 
# print("Model for correlation values without floor effects with a transformation to make data more Normal\n")
# m5 <- lmer(RespCorr_RawTXF ~ CF_kHz + StimSNR + SRgroup + Q10local + TTR_dB + (1|ChinID), data=dat)
# Anova(m5, test.statistic='F')
# # anova(m5)
# plot(m5, main = "Transformed Raw corr")
# plot(m5)
# 
# print("Model for correlation values with floor effects without any transformation\n")
# m4 <- lmer(sigRespCorrRaw ~ sigCF_kHz + sigStimSNR + sigSRgroup + sigQ10local + sigTTR_dB + (1|sigChinID), data=dat)
# Anova(m4, test.statistic='F')
# # anova(m4)
# # plot(m4, main = "Raw corr")
# m4_pred <- fitted(m4)
# # plot(dat$RespCorrRaw, m4_pred)
