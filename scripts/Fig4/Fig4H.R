library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

fils = list.files(path = "../../data/Fig4", pattern = "Fig4H_", full.names = T)
# 
# > fils
# [1] "../../data/Fig4/Fig4H_L428_IL13.csv"     "../../data/Fig4/Fig4H_L428_TARC.csv"    
# [3] "../../data/Fig4/Fig4H_UHO1_IL3_IL13.csv" "../../data/Fig4/Fig4H_UHO1_IL3_TARC.csv"
# [5] "../../data/Fig4/Fig4H_UHO1_IL5_IL13.csv" "../../data/Fig4/Fig4H_UHO1_IL5_TARC.csv"

t_compare_pacritinib = function(data,
                              clone_col = "Clones",
                              value_col = 2,
                              p_adjust_method = "holm") {# should use this default when number of comparisons >= 5
  # Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70. https://www.jstor.org/stable/4615733.
  data[[clone_col]] = as.character(data[[clone_col]])
  clones = unique(data[[clone_col]])
  
  # find pairs *_pacritinib vs base
  pac_pairs = lapply(clones, function(x) {
    if (grepl("_pacritinib$", x)) {
      base = sub("_pacritinib$", "", x)
      if (base %in% clones) return(c(base, x))
    }
    NULL
  })
  pac_pairs = pac_pairs[!vapply(pac_pairs, is.null, logical(1))]
  
  res_list = lapply(pac_pairs, function(pair) {
    base = pair[1]; pac = pair[2]
    y_base = data[[value_col]][data[[clone_col]] == base]
    y_pac  = data[[value_col]][data[[clone_col]] == pac]
    n_base = length(y_base); n_pac = length(y_pac)
    y_all = c(y_base, y_pac)
    total_n = length(y_all)
    
    # --- Classical tests ---
    ttest_p   = t.test(y_pac, y_base, var.equal = FALSE)$p.value

    obs_diff = mean(y_pac) - mean(y_base)
    combs = combn(total_n, n_pac)
    perm_diff = apply(combs, 2, function(idx) mean(y_all[idx]) - mean(y_all[-idx]))
    perm_p = mean(abs(perm_diff) >= abs(obs_diff))
    
    data.frame(
      base = base,
      pacritinib = pac,
      n_base = n_base,
      n_pac = n_pac,
      mean_base = mean(y_base),
      mean_pac = mean(y_pac),
      ttest_p = ttest_p,
      stringsAsFactors = FALSE
    )
  })
  
  res_df = do.call(rbind, res_list)
  
  # Adjust p-values
  res_df$ttest_p_adj  = p.adjust(res_df$ttest_p, method = p_adjust_method)
 
  return(res_df)
}

final_tests = lapply(fils, function(xx){
  dat = read.csv(xx)
  dat[,2] = log2(dat[,2])
  res = compare_pacritinib(dat)
  write.csv(res, gsub(".csv", "_log2scale_Ttest.csv", xx))
})

