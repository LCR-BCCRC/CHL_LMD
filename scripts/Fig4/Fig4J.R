## R code to create figure 4J in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
############## input data ################ 
input_Figure_4J <- read_tsv("../../data/Fig4/FigJ_Affy_CCL17_CSF2RB_MUT.tsv")
##########################
file_name <- "Figure_4J_CCL17_CSF2RBMUT.pdf"
my_plot <- ggplot(Input_Figure_4J, aes(x = CSF2RB_MUT, y = CCL17)) +
  geom_boxplot(alpha = 0.3, color = "black", outlier.shape = NA) +
  geom_beeswarm(aes(fill = CSF2RB_MUT),
                shape = 21, size = 3.5, cex = 2, alpha = 0.6) +
  stat_compare_means() +
  scale_fill_manual(values = c("steelblue", "coral3")) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)) +
  ylab(bquote(italic("CCL17") ~ "transcript level")) +
  xlab(bquote(italic("CSF2RB") ~ "status"))
file_name <- paste(directory_Plots_OI,
                   "CSF2RB_Mut_vs_CCL17_paper.pdf",
                   sep = "/")
print(my_plot)
pdf(file_name)
print(my_plot)
dev.off()
