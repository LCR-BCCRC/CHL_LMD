
#--------------------------------------------------

# Load required packages
library(ggplot2)
library(dplyr)
library(car) 
library(emmeans) 

# Main analysis function
analyze_experiment_data = function(csv_files, 
                                   alpha_normality = 0.05,
                                   alpha_significance = 0.05) {
  
  if (length(csv_files) == 0) {
    stop("No CSV files starting with 'Fig' found in the directory.")
  }
  
  results = list()
  
  for (file in csv_files) {
    cat("Processing:", basename(file), "\n")
    
    tryCatch({
      # Read and prepare data
      data = read.csv(file, stringsAsFactors = FALSE)
      
      # Assume first column is group/factor, second column is values
      colnames(data) = c("Group", "Value")
      data$Group = as.factor(data$Group)
      
      # revised on 20251027
      if("WT" %in% data$Group){
        reference_level = "WT"
      }else{
        # Set reference level to group with shortest name
        reference_level = levels(data$Group)[which.min(nchar(levels(data$Group)))]
      }
      
      data$Group = relevel(data$Group, ref = reference_level)
      
      # Remove missing values
      data = na.omit(data)
      
      # Store file-specific results
      file_result = list(
        filename = basename(file),
        data = data,
        normality_results = list(),
        transformation_used = "raw",  # Default assumption
        test_type = NULL,
        comparisons = list()
      )
      
      # Check normality on raw data
      normality_raw = shapiro.test(data$Value)
      file_result$normality_results$raw = normality_raw
      
      # Determine transformation based on normality test
      if (normality_raw$p.value > alpha_normality) {
        # Raw data is normal - use raw data and stop here
        file_result$transformation_used = "raw"
        data$Value_transformed = data$Value
      } else {
        # Raw data NOT normal - try transformations
        cat("  Raw data not normal, attempting transformations...\n")
        
        # Try log2 transformation (only if all values > 0)
        if (all(data$Value > 0)) {
          data$Value_log2 = log2(data$Value)
          normality_log2 = shapiro.test(data$Value_log2)
          file_result$normality_results$log2 = normality_log2
          
          if (normality_log2$p.value > alpha_normality) {
            file_result$transformation_used = "log2"
            data$Value_transformed = data$Value_log2
            cat("  Using log2 transformation (data now normal)\n")
          } else {
            file_result$transformation_used = "nonparametric"
            data$Value_transformed = data$Value
            cat("  Using nonparametric approach (log2 also not normal)\n")
          }
        } else {
          file_result$transformation_used = "nonparametric"
          data$Value_transformed = data$Value
          cat("  Using nonparametric approach (negative values prevent log transform)\n")
        }
      }
      
      # Create histograms
      plot_list = create_histograms(data, file_result$transformation_used, basename(file))
      file_result$plots = plot_list
      
      # Perform statistical tests
      if (length(unique(data$Group)) == 2) {
        # Two groups case
        file_result = perform_two_group_test(data, file_result, 
                                             reference_level, alpha_significance)
      } else {
        # Multiple groups case
        file_result = perform_multi_group_test(data, file_result, 
                                               reference_level, alpha_significance)
      }
      
      results[[basename(file)]] = file_result
      
    }, error = function(e) {
      cat("Error processing", basename(file), ":", e$message, "\n")
      # Store error information
      results[[basename(file)]] = list(
        filename = basename(file),
        error = e$message
      )
    })
  }
  
  return(results)
}


# Create histogram plots
create_histograms = function(data, transformation, filename) {
  plots = list()
  
  # Raw data histogram
  p1 = ggplot(data, aes(x = Value)) +
    geom_histogram(bins = 10, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = paste("Raw Data -", filename),
         x = "Value", y = "Frequency") +
    theme_minimal()
  
  # Transformed data histogram (if applicable)
  if (transformation == "log2") {
    p2 = ggplot(data, aes(x = Value_log2)) +
      geom_histogram(bins = 10, fill = "lightgreen", color = "black", alpha = 0.7) +
      labs(title = paste("Log2 Transformed -", filename),
           x = "Log2(Value)", y = "Frequency") +
      theme_minimal()
  } else {
    p2 = NULL
  }
  
  return(list(raw = p1, transformed = p2))
}

# Two-group comparison function
perform_two_group_test = function(data, file_result, reference_level, alpha) {
  groups = unique(data$Group)
  
  if (length(groups) != 2) {
    stop("This function is for two-group comparisons only")
  }
  
  # Ensure reference level exists
  if (!reference_level %in% groups) {
    reference_level = groups[1]
    cat("Reference level not found, using", reference_level, "as reference\n")
  }
  
  test_group = setdiff(groups, reference_level)
  
  ref_data = data$Value_transformed[data$Group == reference_level]
  test_data = data$Value_transformed[data$Group == test_group]
  
  if (file_result$transformation_used %in% c("raw", "log2")) {
    # Parametric test (t-test)
    test_result = t.test(test_data, ref_data)
    file_result$test_type = "t-test"
    file_result$comparisons[[paste(test_group, "vs", reference_level)]] = list(
      statistic = test_result$statistic,
      p.value = test_result$p.value,
      significant = test_result$p.value < alpha
    )
  } else {
    # Non-parametric test (Wilcoxon)
    test_result = wilcox.test(test_data, ref_data)
    file_result$test_type = "wilcoxon"
    file_result$comparisons[[paste(test_group, "vs", reference_level)]] = list(
      statistic = test_result$statistic,
      p.value = test_result$p.value,
      significant = test_result$p.value < alpha
    )
  }
  
  return(file_result)
}

perform_multi_group_test = function(data, file_result, reference_level, alpha) {
  groups = unique(data$Group)
  
  if (!reference_level %in% groups) {
    reference_level = groups[1]
    cat("Reference level not found, using", reference_level, "as reference\n")
  }
  
  # Set reference level for the factor
  data$Group = relevel(data$Group, ref = reference_level)
  
  if (file_result$transformation_used %in% c("raw", "log2")) {
    # Parametric: Linear model with normal errors
    file_result$test_type = "linear_model"
    
    # Fit linear model
    lm_model = lm(Value_transformed ~ Group, data = data)
    file_result$model = summary(lm_model)
    file_result$anova = anova(lm_model)
    
    # Overall F-test p-value
    overall_p = anova(lm_model)$`Pr(>F)`[1]
    
    if (overall_p < alpha) {
      # Get coefficients and p-values for each group vs reference
      model_summary = summary(lm_model)
      coefficients = model_summary$coefficients
      
      comparisons = list()
      for (i in 2:nrow(coefficients)) { # Skip intercept (reference level)
        group_name = gsub("Group", "", rownames(coefficients)[i])
        comparison_name = paste(group_name, "vs", reference_level)
        
        comparisons[[comparison_name]] = list(
          estimate = coefficients[i, 1],
          statistic = coefficients[i, 3], # t-value
          p.value = coefficients[i, 4],
          significant = coefficients[i, 4] < alpha
        )
      }
      
      file_result$comparisons = comparisons
    }
    
  } else {
    # Non-parametric: Robust regression or quantile regression
    file_result$test_type = "robust_regression"
    
    # Option 1: Robust regression using M-estimation (rlm from MASS)
    if (require(MASS)) {
      rlm_model = rlm(Value_transformed ~ Group, data = data)
      model_summary = summary(rlm_model)
      
      # For robust regression, we need to approximate p-values
      # Calculate approximate t-tests using the estimated coefficients and SEs
      coefficients = coef(model_summary)
      t_stats = coefficients[, 3]  # t-values
      p_values = 2 * pt(abs(t_stats), df = model_summary$df[2], lower.tail = FALSE)
      
    } else {
      # Option 2: Quantile regression (median regression) using quantreg
      if (require(quantreg)) {
        rq_model = rq(Value_transformed ~ Group, data = data, tau = 0.5)
        model_summary = summary(rq_model, se = "iid") # iid for independent errors
        
        coefficients = coef(model_summary)
        t_stats = coefficients[, 3]  # t-values
        p_values = coefficients[, 4] # p-values
      } else {
        stop("Please install either 'MASS' or 'quantreg' package for non-parametric analysis")
      }
    }
    
    file_result$model = model_summary
    
    # Extract comparisons vs reference
    comparisons = list()
    for (i in 2:length(t_stats)) { # Skip intercept (reference level)
      group_name = gsub("Group", "", names(t_stats)[i])
      comparison_name = paste(group_name, "vs", reference_level)
      
      comparisons[[comparison_name]] = list(
        estimate = coefficients[i, 1],
        statistic = t_stats[i],
        p.value = p_values[i],
        significant = p_values[i] < alpha
      )
    }
    
    file_result$comparisons = comparisons
  }
  
  return(file_result)
}


# Function to generate summary report
generate_summary_report = function(results) {
  cat("EXPERIMENTAL DATA ANALYSIS SUMMARY\n")
  cat("===================================\n\n")
  
  for (file_name in names(results)) {
    result = results[[file_name]]
    
    cat("File:", result$filename, "\n")
    
    if (!is.null(result$error)) {
      cat("ERROR:", result$error, "\n\n")
      next
    }
    
    cat("Transformation used:", result$transformation_used, "\n")
    cat("Test type:", result$test_type, "\n")
    
    # Normality results
    cat("Normality tests:\n")
    for (trans in names(result$normality_results)) {
      p_val = result$normality_results[[trans]]$p.value
      cat("  ", trans, ": p =", round(p_val, 4), 
          ifelse(p_val > 0.05, "(Normal)", "(Non-normal)"), "\n") # this line should be changed as well, from 0.1 to 0.05
    }
    
    # Statistical results
    if (length(result$comparisons) > 0) {
      cat("Comparisons vs", reference_level, ":\n")
      for (comp_name in names(result$comparisons)) {
        comp = result$comparisons[[comp_name]]
        cat("  ", comp_name, ": p =", round(comp$p.value, 4),
            ifelse(comp$significant, "*", ""), "\n")
      }
    } else {
      cat("No significant differences found or no comparisons performed.\n")
    }
    cat("\n" + rep("-", 50) + "\n\n")
  }
}

# Function to save all plots
save_all_plots = function(results, output_directory) {
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  for (file_name in names(results)) {
    result = results[[file_name]]
    
    if (!is.null(result$error) || is.null(result$plots)) {
      next
    }
    # Save raw data histogram
    ggsave(filename = file.path(output_directory, 
                                paste0(tools::file_path_sans_ext(result$filename), 
                                       "_raw_histogram.png")),
           plot = result$plots$raw +
             theme_bw() +  # White background
             theme(panel.background = element_rect(fill = "white"),
                   plot.background = element_rect(fill = "white")),
           width = 8, height = 6, dpi = 300)
    
    # Save transformed data histogram if exists
    if (!is.null(result$plots$transformed)) {
      ggsave(filename = file.path(output_directory, 
                                  paste0(tools::file_path_sans_ext(result$filename), 
                                         "_log2_histogram.png")),
             plot = result$plots$transformed +
               theme_bw() +  # White background
               theme(panel.background = element_rect(fill = "white"),
                     plot.background = element_rect(fill = "white")),
             width = 8, height = 6, dpi = 300)
    }
    
  }
  
  
}

# Test with specific data 
test_specific_file = function(file_path, reference_level = "WT") {
  data = read.csv(file_path)
  colnames(data) = c("Group", "Value")
  data$Group = as.factor(data$Group)
  
  cat("Data structure:\n")
  print(table(data$Group))
  cat("\n")
  
  # Test normality
  normality_raw = shapiro.test(data$Value)
  cat("Raw data normality test p-value:", normality_raw$p.value, "\n")
  
  if (all(data$Value > 0)) {
    data$Value_log2 = log2(data$Value)
    normality_log2 = shapiro.test(data$Value_log2)
    cat("Log2 transformed normality test p-value:", normality_log2$p.value, "\n")
  }
  
  return(data)
}

#----------------------------------------

library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

fils = c("../../data/Fig4/Fig4D.csv", "../../data/Fig4/Fig4G_IL13.csv", "../../data/Fig4/Fig4G_TARC.csv")

results = analyze_experiment_data(fils)

sink("Fig4D_4G_results.txt")
results
sink()
