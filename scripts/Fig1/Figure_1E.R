## R code to create figure 1E in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
######################################################### 
############ load packages ############
######################################################### 
library(readr)
library(dplyr)
library(magrittr)
library(stats)
library(ggplot2)
library(cowplot)
#########################################################
############ Input data ############
#########################################################
input_Figure_1E <- read_tsv("../../data/Fig1/Input_Figure_1E.tsv")
colors_custom = c("#FF0000","#838B8B")
#########################################################
############ Statistics ############
#########################################################
stats_STAT6_MUT <- input_Figure_1E %>%
  dplyr::group_by(AGE_group45) %>%
dplyr::summarise(number_MUT = sum(STAT6_MUT),
          number_WT = n() - number_MUT,
          percentage_MUT = round(100*sum(STAT6_MUT)/n(),digits = 1),
          percentage_WT = 100 - percentage_MUT) 
############ Perform Fisher's exact test ############
matrix_STAT6_MUT <- stats_STAT6_MUT %>%
  dplyr::select(AGE_group45, number_MUT, number_WT) %>%
  tibble::column_to_rownames("AGE_group45") %>%
  as.matrix()
fisher_results <- fisher.test(matrix_STAT6_MUT)
label_for_plot <- paste("p=",round(fisher_results$p.value,4),
                        "\n(Fisher's exact test)",sep = "")
#########################################################
############ Create and save plot ############
#########################################################
stats_STAT6_MUT <- stats_STAT6_MUT[,c("AGE_group45","percentage_WT","percentage_MUT")] %>%
  pivot_longer(-c(AGE_group45),
               names_to = "Subset",
               values_to = "percentage")
stats_STAT6_MUT$Subset <- str_replace(stats_STAT6_MUT$Subset,
                                      "percentage_","")
stats_STAT6_MUT$Subset <- factor(stats_STAT6_MUT$Subset,levels = c("WT","MUT"))
stats_STAT6_MUT %>%
  ggplot() +
  geom_col(aes(Subset,percentage,
               fill= Subset)) +
  facet_grid(~ factor(AGE_group45,levels = c("YOUNG","OLD")))  +
  scale_fill_manual(values = c("MUT" = colors_custom[1],
                               "WT" = colors_custom[2])) +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  geom_text(data = data.frame(Subset = "MUT", percentage = 70, AGE_group45 = "YOUNG", 
                              label = label_for_plot),
              aes(x = Subset, y = percentage, label = label_for_plot),
            size = 4) 

ggsave("STAT6_muts_AGE_plot.pdf",height = 7, width =6)
