
#----------------------libraries--------------------
library(infotheo)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)

#------------------- functions ---------------

mutual_info_importance = function(data, response) {
  
  ## remove UNCLASS on 20241015
  tmp = which(response == "UNCLASS")
  if(length(tmp) > 0){
    response = response[-tmp]
    data = data[-tmp,]
  }
  
  # Ensure response is a factor for mutual information calculation
  if (!is.factor(response)) {
    response = as.factor(response)
  }
  
  levels_response = levels(response)  # Get response levels
  
  num_levels = length(levels_response)
  
  mutual_info_matrix = sapply(names(data), function(feature) {
    feature_data = data[[feature]]
    
    if (is.numeric(feature_data)) {
      # Check if the numeric feature is binary
      if (all(feature_data %in% c(0, 1))) {
        feature_data = as.factor(feature_data)  # Binary features as factors
      } else {
        feature_data = discretize(feature_data)  # Discretize continuous features
      }
    } else {
      feature_data = as.factor(feature_data)  # Convert categorical features to factors
    }
  
    sapply(levels_response, function(level) {
      # Create a binary response for the current level (1 vs. all)
      binary_response = as.factor(response == level)
      mutinformation(feature_data, binary_response, method = "mm") # change on 20241022: set method to mm
    })

  })

  mutual_info_df = as.data.frame(t(mutual_info_matrix))
  colnames(mutual_info_df) = levels_response
  rownames(mutual_info_df) = names(data)
  
  return(mutual_info_df)

}

test_feature_classification_association = function(feature_binary_col, class_col, data) {
  
  # Assume all data are 0 or 1, and set levels as c("1","0")
  data[[feature_binary_col]] = as.factor(data[[feature_binary_col]])
  data[[feature_binary_col]] = factor(data[[feature_binary_col]], levels = c("1", "0"))
  
  # Exclude rows where classification is "UNCLASS"
  data_filtered = data %>% filter(.data[[class_col]] != "UNCLASS")
  
  # Perform chi-squared test for independence
  contingency_table = table(data_filtered[[feature_binary_col]], data_filtered[[class_col]])
  chi_test = chisq.test(contingency_table, simulate.p.value = TRUE)
  
  # add Fisher exact test
  fisher_test = fisher.test(contingency_table)
  
  # Calculate odds ratios for enrichment
  odds_ratios = apply(contingency_table, 2, function(col) {
    # Build a 2x2 table for each class level
    class_yes = col[1] # Number of level1 for this class
    class_no = col[2]  # Number of level2 for this class
    other_yes = sum(contingency_table[1, ]) - class_yes # Other classes for level 1
    other_no = sum(contingency_table[2, ]) - class_no   # Other classes for level 2
    odds_ratio = (class_yes / class_no) / (other_yes / other_no)
    return(odds_ratio)
  })
  
  ## I actually need log odds ratio as well
  log_odds_ratio = log(odds_ratios)
  
  # Return results: chi-squared test p-value and odds ratios
  result = c(
    FisherExactTest_p_value = fisher_test$p.value,
    simulationTest_p_value = chi_test$p.value,
    odds_ratios = unlist(odds_ratios),
    log_odds_ratio
  )
  return(result)
}

mutual_infoPlot_topN_italic = function(mutual_info, titleName, Feature, topN = 5,
                                       classColor = c("#A6FFB6", "#33A1F6", "#D16303","#9370DB")){
  mutual_info_data = cbind(Feature, mutual_info)
  mutual_info_data = data.frame(mutual_info_data)
  mutual_info_data[,2] = as.numeric(mutual_info_data[,2])
  
  # Remove NA
  mutual_info_data = mutual_info_data[!is.na(mutual_info_data[,2]),]
  
  # Order and select top N
  mutual_info_data = mutual_info_data[order(-mutual_info_data$mutual_info), ][1:topN, ]
  
  # Transform Feature names
  mutual_info_data$Feature = gsub("^_", "", mutual_info_data$Feature)                     # Remove leading underscores
  mutual_info_data$Feature = gsub("_", " ", mutual_info_data$Feature)                     # Replace underscores with spaces
  mutual_info_data$Feature = gsub("CNA|MUT", "", mutual_info_data$Feature)               # Remove "CNA" and "MUT"
  
  # Italicize names starting with letters, leave others plain
  mutual_info_data$Feature = ifelse(grepl("^[A-Za-z]", mutual_info_data$Feature),
                                    paste0("italic('", mutual_info_data$Feature, "')"),
                                    mutual_info_data$Feature)
  
  # Lollipop Plot with Flipped Coordinates and Adjusted Labels
  p0 = ggplot(mutual_info_data, aes(x = reorder(Feature, mutual_info), y = mutual_info )) +
    geom_point(color = classColor, size = 3) +
    geom_segment(aes(x = Feature, xend = Feature, y = 0, yend = mutual_info), color = classColor) +
    labs(title = titleName, x = "Features", y = "Mutual Information") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.ticks.length.y = unit(0, "inches")
    ) +
    coord_flip()
  
  # Parse plotmath expressions for axis labels
  p0 = p0 + scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      if (grepl("^italic\\(.*\\)$", label)) {
        parse(text = label)
      } else {
        label
      }
    })
  })
  
  return(list(p0, mutual_info_data))
}


#-----------------------------------------------------------------------------------------------------------------------------------------------------
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

#---------------------------------------------------


#------------------------- LMD cohort --------------------------
classres = read.csv("../../data/Fig6/LMD_predictions.csv", row.names = 1) # nruns = 100 was set for HLGen

## define levels
classres$prediction = factor(classres$prediction, levels = c("CST", "CN913", "STB", "CN2P")) # NO UNCLASS for LMD cohort

## read in data as well
dat = read.csv("../../data/LMD_data.csv", header = T, row.names = 1, stringsAsFactors = F)

classres = classres[rownames(dat),]

## mutual info
minfo = mutual_info_importance(data = dat, response = classres$prediction) 

## I need to combine them together
vars = colnames(dat)[1:dim(dat)[2]]
dat = cbind(dat, classres) 

## test
c_test = t(sapply(vars, test_feature_classification_association, class_col = "prediction", data = dat))
colnames(c_test)[7:10] = paste0("LogOR_", colnames(c_test)[7:10])

## actually, the removing is to set them to NA
info = cbind(c_test, minfo[rownames(c_test),])

colnames(info)[11:14] = paste0("Mutual_info_",colnames(info)[11:14]) 

# > colnames(info)
# [1] "FisherExactTest_p_value" "simulationTest_p_value"  "odds_ratios.CST"         "odds_ratios.CN913"       "odds_ratios.STB"        
# [6] "odds_ratios.CN2P"        "LogOR_CST"               "LogOR_CN913"             "LogOR_STB"               "LogOR_CN2P"             
# [11] "Mutual_info_CST"         "Mutual_info_CN913"       "Mutual_info_STB"         "Mutual_info_CN2P"

info$Mutual_info_CST = ifelse(info$odds_ratios.CST < 1, NA, info$Mutual_info_CST)
info$Mutual_info_CN913 = ifelse(info$odds_ratios.CN913 < 1, NA, info$Mutual_info_CN913)
info$Mutual_info_STB = ifelse(info$odds_ratios.STB < 1, NA, info$Mutual_info_STB)
info$Mutual_info_CN2P = ifelse(info$odds_ratios.CN2P < 1, NA, info$Mutual_info_CN2P)

titleNames = c("CST", "CN913", "STB", "CN2P")
allplots = list()
impscores = list()
colors = c("#A6FFB6", "#33A1F6", "#D16303","#9370DB")

for (i in 1:4){
  tmp = mutual_infoPlot_topN_italic(info[,i+10], titleNames[i], Feature = rownames(info), topN = 5, classColor = colors[i])
  allplots[[i]] = tmp[[1]]
  impscores[[i]] = tmp[[2]]
}

combined_plot_cowplot = plot_grid(plotlist = allplots, ncol = 4, nrow = 1)
print(combined_plot_cowplot)
ggsave("Fig6B_LMD_top5positive_mutualInfoPlots.pdf", height = 3, width = 13)



#---------------------------------------------------------------------------
#------------------------- Alig cohort --------------------------
classres = read.csv("../../data/Fig6/Alig_predictions.csv", row.names = 1) # nruns = 100 was set for HLGen

#------------------------------------------------------------------------------------------------
## define levels, this is the most important code line in this file
classres$prediction = factor(classres$prediction, levels = c("CST", "CN913", "STB", "CN2P", "UNCLASS"))
#-----------------------------------------------------------------------------------------------

## read in data as well
dat = read.csv("../../data/Alig_data.csv", header = T, row.names = 1, stringsAsFactors = F)

classres = classres[rownames(dat),]

## mutual info
minfo = mutual_info_importance(data = dat, response = classres$prediction)

## I need to combine them together
vars = colnames(dat)[1:dim(dat)[2]]
dat = cbind(dat, classres) 

#------------------ last step to deal with UNCLASS------------------------------------
# > head(minfo)
#              CST       CN913          STB        CN2P      UNCLASS
# ACTB_MUT   0.002972102 0.013330459 0.000000e+00 0.004419677       0
# ARID1A_MUT 0.000000000 0.006381467 3.824708e-05 0.000000000       0
# B2M_MUT    0.039516820 0.037486942 1.732765e-01 0.002186586       0
# BCL2_MUT   0.000112410 0.006102210 0.000000e+00 0.025891335       0
# BCL7A_MUT  0.053241512 0.088225689 0.000000e+00 0.000000000       0
# BTG1_MUT   0.006463927 0.008239052 0.000000e+00 0.000000000       0
minfo = minfo[, 1:4]

dat$prediction = as.character(dat$prediction)
dat$prediction = ifelse(dat$prediction == "UNCLASS", NA, dat$prediction)
# > table(dat$prediction, useNA = "ifany")
# 
# CN2P CN913   CST   STB  <NA> 
#   16    45   114    66    52 
dat$prediction = factor(dat$prediction, levels = c("CST", "CN913", "STB", "CN2P"))
#------------------------------------------------------

## test
c_test = t(sapply(vars, test_feature_classification_association, class_col = "prediction", data = dat))
colnames(c_test)[7:10] = paste0("LogOR_", colnames(c_test)[7:10])

info = cbind(c_test, minfo[rownames(c_test),])

colnames(info)[11:14] = paste0("Mutual_info_",colnames(info)[11:14]) 

# > colnames(info)
# [1] "FisherExactTest_p_value" "simulationTest_p_value"  "odds_ratios.CST"         "odds_ratios.CN913"       "odds_ratios.STB"        
# [6] "odds_ratios.CN2P"        "LogOR_CST"               "LogOR_CN913"             "LogOR_STB"               "LogOR_CN2P"             
# [11] "Mutual_info_CST"         "Mutual_info_CN913"       "Mutual_info_STB"         "Mutual_info_CN2P"

info$Mutual_info_CST = ifelse(info$odds_ratios.CST < 1, NA, info$Mutual_info_CST)
info$Mutual_info_CN913 = ifelse(info$odds_ratios.CN913 < 1, NA, info$Mutual_info_CN913)
info$Mutual_info_STB = ifelse(info$odds_ratios.STB < 1, NA, info$Mutual_info_STB)
info$Mutual_info_CN2P = ifelse(info$odds_ratios.CN2P < 1, NA, info$Mutual_info_CN2P)

titleNames = c("CST", "CN913", "STB", "CN2P")
allplots = list()
impscores = list()
colors = c("#A6FFB6", "#33A1F6", "#D16303","#9370DB")

for (i in 1:4){
  tmp = mutual_infoPlot_topN_italic(info[,i+10], titleNames[i], Feature = rownames(info), topN = 5, classColor = colors[i])
  allplots[[i]] = tmp[[1]]
  impscores[[i]] = tmp[[2]]
}

combined_plot_cowplot = plot_grid(plotlist = allplots, ncol = 4, nrow = 1)
print(combined_plot_cowplot)
ggsave("Fig6C_Alig_top5positive_mutualInfoPlots.pdf", height = 3, width = 13)


#----------------------------------------------------
