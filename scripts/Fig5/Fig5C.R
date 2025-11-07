
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# read in the updated data with many 0s
dat = read.csv("../../data/Fig5/Fig5C.csv")

# > dat
#   Clones Values
# 1      WT   6.82
# 2      WT   1.69
# 3      WT   0.00
# 4      WT   0.00
# 5      WT   0.00
# 6      WT   8.01
# 7    E788   0.00
# 8    E788  15.60
# 9    E788  22.90
# 10   E788  13.60
# 11   E788   9.92

dat$isNonZero = ifelse(dat$Values > 0.00, 1, 0)

#=========================================================
# normality test

sdat = subset(dat, dat$isNonZero == 1)

# > sdat 
# Clones Values isNonZero
# 1      WT   6.82         1
# 2      WT   1.69         1
# 6      WT   8.01         1
# 8    E788  15.60         1
# 9    E788  22.90         1
# 10   E788  13.60         1
# 11   E788   9.92         1

ntest = shapiro.test(sdat$Values)
# 
# > ntest
# 
# Shapiro-Wilk normality test
# 
# data:  sdat$Values
# W = 0.97987, p-value = 0.9589
# good, again, normality is hold!

library(censReg)

# haha, I need set the factor level first
dat$Clones = factor(dat$Clones, levels = c("WT", "E788"))

tobit_model = censReg(Values ~ Clones, data = dat)
summary(tobit_model)
# 
# > summary(tobit_model)
# 
# Call:
#   censReg(formula = Values ~ Clones, data = dat)
# 
# Observations:
#   Total  Left-censored     Uncensored Right-censored 
# 11              4              7              0 
# 
# Coefficients:
#   Estimate Std. error t value  Pr(> t)    
# (Intercept)  -0.5398     3.8036  -0.142   0.8871    
# ClonesE788   12.2321     5.2017   2.352   0.0187 *  
#   logSigma      2.0799     0.2897   7.180 6.96e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Newton-Raphson maximisation, 7 iterations
# Return code 1: gradient close to zero (gradtol)
# Log-likelihood: -27.72963 on 3 Df

#----------------------------------------------
# In a Tobit regression model accounting for censoring at zero, the ClonesE788 group 
# showed significantly higher values than the reference group (β = 12.2321, SE = 5.2017, t = 2.353, p = 0.0187)
#---------------------------------------------

