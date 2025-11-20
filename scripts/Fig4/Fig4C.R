## R code to create figure 4C in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
#############################################
############## load packages ################
#############################################
library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
#############################################
############## input data ################ 
#############################################
input_Figure_4C <- read_tsv("../../data/Fig3/Input_Figure_4C")
#############################################
############## custom function ##################
#############################################
create_VAF_plot_CSF2RB <- function(input,
                                           directory_PDF,
                                           name_PDF,
                                           title = NULL,
                                           sample = "res_id",
                                           VAF_tumor = "VAF_tumor",
                                           aes_diff = "GENE %in% c('CSF2RB')",
                                           legend = FALSE,
                                           y_max = 100,
                                           fontsize_x_axis = 12,
                                           fontsize_y_axis = 12,
                                           font_x_axis = "plain",
                                           ylab = "Variant allele frequency",
                                           xlab = "",
                                           log_scale_y = "yes",
                                           fontsize_title_x = 12,
                                           fontsize_title_y = fontsize_title_x,
                                           stroke = .4){
My_plot <- ggplot(input,aes_string(sample,VAF_tumor)) +
  geom_beeswarm(aes_string(fill = aes_diff,size = aes_diff,alpha = aes_diff),
                shape = 21,color = "darkgrey",cex = .4,stroke= stroke) +
  geom_boxplot(alpha = 0,lwd = .3) +
theme_cowplot() +
    guides(fill=guide_legend(title="Mutation"),size = "none",alpha = "none") + #only show legend for fill
  scale_fill_manual(values = c("grey","red"),labels = c("Other",expression(italic("CSF2RB")))) +
  scale_size_manual(values = c(2.5,4)) +
  scale_alpha_manual(values = c(0.4,0.9)) +
  theme(axis.text.x = element_text(face =  font_x_axis,size = fontsize_x_axis,angle = 35,hjust=1),
        axis.text.y = element_text(size = fontsize_y_axis),
        axis.title.x = element_text(size = fontsize_title_x),
        axis.title.y = element_text(size = fontsize_title_y),
        legend.title = element_blank(),## removes legend title 
        legend.text.align = 0,
        strip.text = element_text(size = 12),
        legend.position = c(rel(0.85),rel(0.85)),
        legend.background = element_rect(
      colour = "grey", # Color of the border
      fill = "white", # Background fill color of the legend box
      alpha = 0,
      linetype = "dotted", # Line type of the border (e.g., "solid", "dashed", "dotted")
      linewidth = .1 # Thickness of the border line
    ),
    legend.margin = margin(t = 5, r = 10, b = 5, l = 10))  + ## otherwise the line is to tight around the legendtext
  ylab(ylab) +
  xlab(xlab) +
  geom_hline(yintercept = c(0,.1,.25,.50,1),linetype = "dotted",color = "grey",alpha = 0.6) 
if(legend == FALSE){
  My_plot <- My_plot +
    theme(legend.position = "none") 
}

My_plot <- addSmallLegend(My_plot,textSize = 12,pointSize = 1, spaceLegend = 0.3)


if(log_scale_y == "yes"){
  breaks=c(.0001,.001,.01,.1,.25,.5,1)
  My_plot <- My_plot + 
    scale_y_log10(breaks = breaks)
}
if(!is.null(title)){
  My_plot <- My_plot +
    ggtitle(title)
}
name_PDF <- paste(directory_PDF,name_PDF,sep = "")
print(My_plot)
ggsave(name_PDF, width = 12, height = 4)
}
#############################################
############## Change order input data ##################
#############################################
# order by median variant allele frequency (VAF)
order_plot <- Input_Figure_4C %>%
  arrange(desc(median_VAF), res_id)
order_of_samples <- unique(order_plot$res_id)
Input_Figure_4C$res_id <- factor(Input_Figure_4C$res_id,levels = order_of_samples)
#############################################
############### create and save plot #################
#############################################
create_VAF_plot_CSF2RB(Input_Figure_34,
                               directory_PDF = "",
                              legend =  TRUE,
                               name_PDF = "Figure_4C_VAF_distributions_CSF2RB.pdf",
                              fontsize_title_x = 10,
                               fontsize_x_axis = 8,
                               fontsize_y_axis = 8) 
#############################################
