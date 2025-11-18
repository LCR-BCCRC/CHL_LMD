## R code to create figure 1F in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
######################################################### 
############ load packages ############
######################################################### 
library(survival)
library(survminer)
library(GGally)
library(readr)
######################################################### 
############# read input data #############
######################################################### 
Input_figure_1F <- read_tsv("../../data/Fig1/Input_Figure1.tsv")
######################################################### 
############ create custom function ############
######################################################### 
Analyze_survival_PDF <- function(Samples_to_be_analyzed,
         Variables_stratification,
         title_graph = NULL,
         Survival_Times,
         save_pdf = "no",
         save_tiff = "no",  # NEW FLAG
         settings_y_axis_vjust = 0,
         settings_y_axis_margin_l = 0,
         settings_y_axis_margin_r = 0,
         label = NULL,
         colors = NULL,
         directory = NULL,
         x_axis_max = NULL,
         legend.title = " ",
         break.time.by = 2,
         font.main = 20,
         font.submain = 18,
         font.x = 18,
         font.y = 18,
         font.caption = 16,
         font.tickslab = 14,
         legend = c(.8, .3),
         font.legend = c(12, "bold", "black"),
         risk.table.y.text.col = TRUE,
         risk.table.height = 0.3,
         risk.table.fontsize = 4.5,
         risk.table.font = "bold",
         risk.table.pos = "out",
         risk_table_title_fontsize = 12,
         table_axis.text.fontsize = 12,
         surv.median.line = "hv",
         censor.size = 10,
         base_size = 12,
         pval = TRUE,
         pval.size = 5,
         size = 1.5) {
  
  for(j in Variables_stratification) {
    tryCatch({
      for(i in Survival_Times) {
        tryCatch({
          Time <- paste(i, label, sep = "_")
          CODE <- paste("CODE", Time, sep = "_")
          
          if(!is.null(x_axis_max)) {
            xlim <- c(0, x_axis_max)
          } else {
            xlim <- c(0, max(Samples_to_be_analyzed[[Time]], na.rm = TRUE))
          }
          
          f <- paste("Surv(", Time, ",", CODE, ") ~ ", j, sep = "")
          fit <- do.call("survfit", list(as.formula(f), data = as.name(Samples_to_be_analyzed)))
          print(fit)
          
          if(is.null(fit$strata)) {
            labels_ <- "all"
          } else {
            labels_ <- str_replace_all(unique(names(fit$strata)), ".+[//=]", "")
          }
          
          if(is.null(colors)) {
            colors_ <- pal_jco()(length(labels_))
          } else {
            colors_ <- colors
          }
          
          if(is.null(title_graph)) {
            title_graph <- ""
          }
          
          ggsurv <- ggsurvplot(fit,
                               ggtheme = theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(),
                                               panel.background = element_blank(),
                                               axis.title.y = element_text(vjust = settings_y_axis_vjust,
                                                                           size = 8,
                                                                           margin = margin(l = settings_y_axis_margin_l,
                                                                                           r = settings_y_axis_margin_r)),
                                               axis.line = element_line(colour = "black"),
                                               legend.background = element_rect(fill = "transparent"),
                                               legend.box.background = element_blank(),
                                               legend.key = element_blank()),
                               title = title_graph,
                               pval = pval,
                               pval.size = pval.size,
                               size = size,
                               break.time.by = break.time.by,
                               base_size = base_size,
                               base_family = "",
                               font.main = c(font.main, "bold", "black"),
                               font.submain = c(font.submain, "plain", "black"),
                               font.x = c(font.x, "bold", "black"),
                               font.y = c(font.y, "bold", "black"),
                               font.caption = c(font.caption, "plain", "black"),
                               font.tickslab = c(font.tickslab, "plain", "black"),
                               xlim = xlim,
                               legend = legend,
                               font.legend = font.legend,
                               risk.table = TRUE,
                               risk.table.y.text.col = risk.table.y.text.col,
                               risk.table.height = risk.table.height,
                               risk.table.fontsize = risk.table.fontsize,
                               risk.table.font = risk.table.font,
                               risk.table.pos = risk.table.pos,
                               censor.size = censor.size,
                               legend.labs = labels_,
                               tables.theme = theme(plot.title = element_text(size = risk_table_title_fontsize),
                                                    axis.text.y = element_text(size = table_axis.text.fontsize),
                                                    axis.text.x = element_blank(),
                                                    axis.title = element_blank()),
                               surv.median.line = surv.median.line,
                               xlab = "Time (years)",
                               ylab = paste("Survival Probability (", i, ")", sep = ""),
                               palette = colors_,
                               legend.title = legend.title)
          
          # Save plots (PDF, TIFF) if requested
          if(save_pdf == "yes" || save_tiff == "yes") {
            if(is.null(directory)) {
              file_base <- paste(j, "_", i, sep = "")
            } else {
              file_base <- paste(directory, j, "_", i, sep = "")
            }
            
            # Save PDF
            if(save_pdf == "yes") {
              pdf(paste0(file_base, ".pdf"), width = 10, height = 8, onefile = FALSE)
              print(ggsurv)
              dev.off()
            }
            
            # Save TIFF
            if(save_tiff == "yes") {
              tiff(paste0(file_base, ".tiff"), width = 10, height = 8, units = "in", res = 300, compression = "lzw")
              print(ggsurv)
              dev.off()
            }
          }
          
          print(ggsurv)
          
        }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
      }
    }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
  }
}
#########################################################  
############input parameters for survival ############
#########################################################  
Variables_stratification <- c("STAT6_mut_gain")
Survival_times <- c("PFS")
input_survival <- "Input_figure_1F"
colors_custom = c("#FF0000","#838B8B")

#########################################################  
############ plot and save figure ############
#########################################################
Analyze_survival_PDF(Samples_to_be_analyzed = input_survival,
                 Variables_stratification = Variables_stratification,
                 size =2,
                 save_pdf = "yes",
                 save_tiff = "yes",
                 colors = colors_custom,
                 Survival_Times = Survival_Times,
                  x_axis_max = 10,
                 legend = c(0.8, 0.1),
                             font.legend = c(9, "bold", "black"),
                 label = "Survival"
                )
