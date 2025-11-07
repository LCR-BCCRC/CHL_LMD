
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# read in data with many 0s

dat = read.csv("../../data/Fig4/Fig4K_CSF2RB_mut_TARC.csv")
sdat = subset(dat, dat$CCL17_percentage_HL_cohort > 0)
shapiro.test(as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-WT")$CCL17_percentage_HL_cohort))
shapiro.test(as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-Mutant")$CCL17_percentage_HL_cohort))

# > shapiro.test(as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-WT")$CCL17_percentage_HL_cohort))
# 
# Shapiro-Wilk normality test
# 
# data:  as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-WT")$CCL17_percentage_HL_cohort)
# W = 0.74448, p-value = 1.248e-09
# 
# > shapiro.test(as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-Mutant")$CCL17_percentage_HL_cohort))
# 
# Shapiro-Wilk normality test
# 
# data:  as.vector(subset(sdat, sdat$CSF2RB_MUT == "CSF2RB-Mutant")$CCL17_percentage_HL_cohort)
# W = 0.7782, p-value = 0.0002325
## both of them are sig, did not pass normality check

dat$isNonZero = ifelse(dat$CCL17_percentage_HL_cohort > 0.00, 1, 0)

# > table(dat$CSF2RB_MUT, dat$isNonZero)
# 
#                0  1
# CSF2RB-Mutant  0 22
# CSF2RB-WT     13 69

# Unconditional Test: Barnard's test       
library(Exact)
# install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')
exact.test(data = as.matrix(table(dat$CSF2RB_MUT, dat$isNonZero)), method = "csm") # this is true Barnard's test 
# > exact.test(data = as.matrix(table(dat$CSF2RB_MUT, dat$isNonZero)), method = "csm") # this is true Barnard's test 
# 
# CSM Exact Test
# 
# data:  0 out of 22 vs. 13 out of 82
# test statistic = NA, first sample size = 22, second sample size = 82, p-value = 0.04294
# alternative hypothesis: true difference in proportion is not equal to 0
# sample estimates:
#   difference in proportion 
# -0.1585366 

# Barnard's CSM test found a statistically significant difference in proportions (p = 0.043) between the two groups, 
# with the mutation group showing 0/22 zero percentage compared to 13/82 zero percentage rate in the WT group, 
# representing a 15.9% difference in proportions with 0 percentage.
# 
# ----------------------------
# References:
# 
# Barnard, G.A. (1945) A new test for 2x2 tables. Nature, 156, 177
# Barnard, G.A. (1947) Significance tests for 2x2 tables. Biometrika, 34, 123–138

