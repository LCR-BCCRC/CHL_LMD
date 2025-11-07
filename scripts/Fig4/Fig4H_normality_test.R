library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

fils = list.files(path = "../../data/Fig4", pattern = "Fig4H_", full.names = T)
# 
# > fils
# [1] "../../data/Fig4/Fig4H_L428_IL13.csv"     "../../data/Fig4/Fig4H_L428_TARC.csv"    
# [3] "../../data/Fig4/Fig4H_UHO1_IL3_IL13.csv" "../../data/Fig4/Fig4H_UHO1_IL3_TARC.csv"
# [5] "../../data/Fig4/Fig4H_UHO1_IL5_IL13.csv" "../../data/Fig4/Fig4H_UHO1_IL5_TARC.csv"

#===============================================
## Overall normality test and export results
#===============================================

sink("normalityTest_Fig4H.txt")

# It's not a good idea to use lapply here; use a for loop instead
for (i in 1:length(fils)){
  xdat = read.csv(fils[i])
  cat("\n")
  cat("=====================================\n\n")
  cat(fils[i])
  print(shapiro.test(xdat[,2]))
  cat("\n")
}

sink()
# Only one dataset is marginal; all others clearly failed the test
# ../data/Fig4/Fig4H_L428_TARC.csv
# Shapiro-Wilk normality test
# 
# data:  xdat[, 2]
# W = 0.93286, p-value = 0.05852


# Try again after log2 transformation

sink("normalityTest_log2transformed_Fig4H.txt")

# It's not a good idea to use lapply here; use a for loop instead
for (i in 1:length(fils)){
  xdat = read.csv(fils[i])
  cat("\n")
  cat("=====================================\n\n")
  cat(fils[i])
  print(shapiro.test(log2(xdat[,2])))
  cat("\n")
}

sink()

# The same dataset that was marginal in raw data now passes
# ../data/Fig4/Fig4H_L428_TARC.csv
# Shapiro-Wilk normality test
# 
# data:  log2(xdat[, 2])
# W = 0.94179, p-value = 0.1017

# Another dataset also passed

# =====================================
#   
# ../data/Fig4/Fig4H_UHO1_IL3_IL13.csv
# Shapiro-Wilk normality test
# 
# data:  log2(xdat[, 2])
# W = 0.91755, p-value = 0.117

## Therefore, these two datasets are suitable for t-tests using log2-transformed data



#===============================================
## Test normality within each group and export results
#===============================================

sink("withinGroup_normalityTest_Fig4H.txt")

# It's not a good idea to use lapply here; use a for loop instead
for (i in 1:length(fils)){
  xdat = read.csv(fils[i])
  cat("\n")
  cat("=====================================\n\n")
  cat(fils[i])
  print(by(xdat[,2], xdat[,1], shapiro.test))
  cat("\n")
}

sink()

# Try again after log2 transformation

sink("withinGroup_normalityTest_log2transformed_Fig4H.txt")

# It's not a good idea to use lapply here; use a for loop instead
for (i in 1:length(fils)){
  xdat = read.csv(fils[i])
  cat("\n")
  cat("=====================================\n\n")
  cat(fils[i])
  print(by(log2(xdat[,2]), xdat[,1], shapiro.test))
  cat("\n")
}

sink()

# While all other tables passed normality test when we work on each category separately, 
# some group levels from ../data/Fig4/Fig4H_L428_TARC.csv showed deviations from normality, 
# their small sample sizes likely contribute to this variability. 
# Because the overall distributions were acceptable after log2  transformation for this table,
# and to maintain consistency across all comparisons, log2-transformed data were used for all subsequent tests in Fig. 4H

