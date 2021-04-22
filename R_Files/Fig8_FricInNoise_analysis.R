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
dat <- read.table('table_fric_in_noise_withNan.txt', header=TRUE)
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

## Final
print("Final")
m4 <- lm(RespCorrRaw  ~ CF_kHz_log + StimSNR + SRgroup*HearingStatus, data=dat)
Anova(m4, test.statistic='F')

## deq check 
m1_fricIN <- lmer(RespCorrRaw ~ CF_kHz_log + StimSNR + Threshold_dBSPL + Q10local + TTR_dB + (1|ChinID), data=dat)
anOut_f1_hf <- Anova(m1_fricIN, test.statistic='F')
print(anOut_f1_hf)
var_out <- get_deq(anOut_f1_hf)
print(var_out)