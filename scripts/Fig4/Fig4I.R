## R code to create figure 4I in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
############## input data ################ 
input_Figure_4I <- read_tsv("../../data/Fig4/Fig4I_cHL_tumors_TARC_CSF2RB_MUT.tsv")
##########################
input_Figure_4I %>%
  ggplot(aes(x = CCL17_percentage_HL_cohort, fill = CSF2RB_MUT)) +
  geom_bar(aes(y = after_stat(prop), group = CSF2RB_MUT), position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_cowplot() +
  scale_fill_manual(values = c("steelblue2","coral3"),
                    labels = c("CSF2RB-WT", "CSF2RB-mutated")) +
  xlab("TARC positivity (percentage of cells)") +
  ylab("proportion of cases") +
  theme(legend.title = element_blank()) +
  ggtitle("HL-LMD cohort")

ggsave("CSF2RB_MUT_CCL17_perc_distr.pdf")
#############
