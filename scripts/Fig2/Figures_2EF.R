## R code to create figures 2.E and 2.F in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readr)
############## define colors ################ 
colors_HLGen <- c("CST" = "#A6FFB6",
                  "CN913" = "#33A1F6",
                  "STB" = "#D16303",
                  "CN2P" = "#9370DB",
                  "UNCLASS" = "#D3D3D3")


############## input data ################ 
input_Figures_2EF_ssgsea <- read_tsv("../../data/Fig2/input_Figures_2EF_ssgsea.tsv")
top_signatures_ssgsea <- read_tsv("../../data/Fig2/input_Figures_2EF_top_signatures.tsv")
##############################
### Set correct order HLGen categories
input_Figures_2EF_ssgsea$HLGen <- factor(input_Figures_2EF_ssgsea$HLGen,c("CST","CN913","STB","CN2P"))
##############################

for(j in c("Figure_2E","Figure_2F")) {
  if(i == "Figure_2E"){signature_OI = "REACTOME_INTERFERON_SIGNALING_ssgsea_lymphoma"}
  if(i == "Figure_2F"){signature_OI = "PEREZ_TP53_TARGETS_ssgsea_lymphoma"}
  ###
  cluster_highlight <- top_signatures_ssgsea$Feature[top_signatures$signature ==signature_OI]
  cluster_highlight <- str_replace(cluster_highlight,"HLGen_binary_k_","")
  #####
  ## define pairs and function for statistical comparison
  #####
  if(cluster_highlight == "CST"){
    pairs_to_compare <- list(c("CST","CN913"),
                             c("CST","STB"),
                             c("CST","CN2P"))}
  if(cluster_highlight == "CN913"){
    pairs_to_compare <- list(c("CST","CN913"),
                             c("CN913","STB"),
                             c("CN913","CN2P"))}
  if(cluster_highlight == "STB"){
    pairs_to_compare <- list(c("CST","STB"),
                             c("CN913","STB"),
                             c("STB","CN2P"))}
  if(cluster_highlight == "CN2P"){
    pairs_to_compare <- list(c("CST","CN2P"),
                             c("CN913","CN2P"),
                             c("STB","CN2P"))}
  #####
  function_stat_compare_means <-    stat_compare_means(comparisons = pairs_to_compare,
                                                       bracket.size = .1,
                                                       size =4,
                                                       step.increase = 0.1,
                                                       method = "wilcox.test",
                                                       tip.length=0.05,
                                                       vjust = 0.3,
                                                       label = "p.signif", 
                                                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                          symbols = c("****", "***", "**", "*", "ns")))
  ###################################
  plot_OI <- ggplot(input_Figures_2EF_ssgsea,aes_string("HLGen",signature_OI)) +
    geom_boxplot(alpha = 0) +
    geom_beeswarm(size = 4.5, alpha = .9, shape =21, cex = 2.2,
                  aes_string(fill = "HLGen",
                             color = "ifelse(HLGen == cluster_highlight,'red','blue')")) +
    function_stat_compare_means +
    ylab("single sample GSEA score") +
    xlab("HLGen") +
    scale_fill_manual(values = colors_HLGen) + 
    scale_color_manual(values = c("white","black"))  +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle(str_replace(g,"_SigDB_Staudt|_ssgsea_lymphoma","")) 
  
  name_pdf <- paste(j,".pdf",sep = "")
  print(plot_OI)  
  ggsave(name_pdf,height = 5, width = 8)
}
