## R code to create figure 4F in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
############## input data ################ 
Input_Fig4F_19_genes_scaled <- read_tsv("../../data/Fig4/Fig4F_L428_19genes_scaled.tsv")
Input_Fig4F_19_genes_scaled <- Input_Fig4F_19_genes_scaled %>%
  column_to_rownames("Gene")
############################
############## Annotatoin for heatmap ##############
df <- data.frame(genotype = str_extract(colnames(Input_Fig4F_19_genes_scaled),"[^-]+"))
df[] <- df[] %>%
  lapply(as.character)
ha_column <- HeatmapAnnotation(df = df,
                               annotation_name_gp = gpar(fontsize = 8),
                               show_legend = TRUE,
                               col = list(genotype = c("Emp" = "grey",
                                                             "WT" = "steelblue3",
                                                       "skip13" = "coral2",
                                                             "E788" = "firebrick"),
                                          treatment = c("none" = "grey",
                                                        "60min" = "coral4")),
                               na_col ="white"
                               )
##############
############## Create heatmap ##############
Heatmap <- Heatmap(Input_Fig4F_19_genes_scaled,
                    name = "z-score",
                   column_title = "L428",
                   column_title_gp = gpar(8),
                   row_dend_side = "right",
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 6, 
                                       fontface = "italic"),
                   top_annotation = ha_column,
                   height = unit(10,"cm"))
pdf("Figure4F_Heatmap_L428_19Genes.pdf",
    width = 4,
    height = 6)
Heatmap <- ComplexHeatmap::draw(Heatmap,show_heatmap_legend = TRUE)
dev.off()
