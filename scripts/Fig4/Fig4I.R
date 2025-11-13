## R code to create figure 4I in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
############## input data ################ 
input_Figure_4I <- read_tsv("../../data/Fig4/Fig4I_DGE_Affy_CSF2RB_MUT.tsv")
##########################
genes_19_signature <- c("KCNJ2", 
"CCL22", 
"SSTR2", 
"TNFRSF9", 
"PM20D2", 
"SLAMF1", 
"CCL17", 
"LTBP2", 
"ATP8B4", 
"SLC26A4", 
"CCND1", 
"RHOV", 
"PTP4A3", 
"IL13", 
"SLCO4A1", 
"KCNA3", 
"NCALD",
"LTA", 
"ZBTB32")
##########################
PLOT_DGEA <- Input_Figure_4I %>%
        ggplot(aes(logFC, -log10(P.Value)))  +
        geom_point(shape = 21, color= "grey", 
                   aes(fill =Gene %in% genes_19_signature,
                       alpha = Gene %in% genes_19_signature,
                       size= Gene %in% genes_19_signature))  +
        scale_fill_manual(values = c("grey","coral4")) +
        scale_alpha_manual(values = c(0.01,0.8)) +
        scale_size_manual(values = c(2,3.5)) +
        geom_vline(xintercept = 0, color = "steelblue2", linetype = "dashed") +
        geom_text_repel(aes(logFC,-log10(P.Value),
                            label=ifelse(Gene %in% genes_19_signature &
                                           -log10(P.Value) > .4,as.character(Gene),"")),
                        fontface = "italic", 
                        size = 4.5,
                        force = 10,
                        segment.alpha = 0.8,
                        segment.color = "grey",
                        segment.size = .3,
                        max.overlaps = Inf) +
  theme_cowplot() +
  theme(legend.position = "none")

 pdf_name <- "Figure_4I_volcano_Affy_CSF2RBMUT.pdf"
print(PLOT_DGEA)
ggsave(pdf_name,height =6, width = 5)
#############
