## R code to create figures 4.A and 4.B in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
################################################
############## load packages ################ 
################################################
library(GAMBLR)
library(tidyverse)
library(readr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
################################################
source("get_gambl_colours.R")
source("create_onco_matrix.R")
source("prettyOncoplot.R")
################################################
############## Input files ###############
################################################
lmd_maf <- read_tsv("../../data/Fig2/Input_Figure_2A.maf")

# Tally of mutation types per gene
mut_gt_5 <- lmd_maf %>%
    distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    distinct() %>%
    count(Hugo_Symbol) %>%
    filter(n >= 10) %>%
    pull(Hugo_Symbol)

mut_prop <- lmd_maf %>%
    filter(Hugo_Symbol %in% mut_gt_5) %>%
    mutate(Variant_Classification = case_when(
        str_detect(Variant_Classification, "Nonsense|Frame_Shift") ~ "Truncating",
        Variant_Classification == "Splice_Site" ~ "Splice Site",
        Variant_Classification == "Missense_Mutation" ~ "Missense",
        str_detect(Variant_Classification, "In_Frame") ~ "In Frame In/Del",
        Variant_Classification %in% c("Translation_Start_Site", "Nonstop_Mutation") ~ str_replace_all(Variant_Classification, "_", " "),
    )) %>%
    filter(!is.na(Variant_Classification)) %>%
    group_by(Hugo_Symbol) %>%
    count(Variant_Classification) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    drop_na()

mut_order <- mut_prop %>%
    select(-n) %>%
    pivot_wider(
        names_from = Variant_Classification,
        values_from = prop
    ) %>%
    arrange(desc(Truncating), desc(Missense)) %>%
    pull(Hugo_Symbol)

mut_prop$Hugo_Symbol <- factor(mut_prop$Hugo_Symbol, levels = rev(mut_order))

ggplot(
    mut_prop,
    aes(
        x = Hugo_Symbol,
        y = prop,
        fill = Variant_Classification
    )
) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
        values = c(
            "Truncating" = unname(get_gambl_colours("mutation")["Nonsense_Mutation"]),
            "Missense" = unname(get_gambl_colours("mutation")["Missense_Mutation"]),
            "Splice Site" = unname(get_gambl_colours("mutation")["Splice_Site"]),
            "In Frame In/Del" = unname(get_gambl_colours("mutation")["In_Frame_Ins"]),
            "Translation Start Site" = unname(get_gambl_colours("mutation")["Translation_Start_Site"]),
            "Nonstop Mutation" = unname(get_gambl_colours("mutation")["Nonstop_Mutation"])
        ),
        name = "Variant Type"
    ) +
    xlab("") +
    ylab("Proportion of Mutations") +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_text(face = "italic"))

ggsave("Figure_4A_mutation_types.pdf", width = 6, height = 6)

lmd_maf$HGVSp_Short <- lmd_maf$Protein_Change
source("pretty_lollipop_plot.R")
pretty_lollipop_plot(
    maf = lmd_maf,
    gene = "CSF2RB"
) +
    cowplot::theme_cowplot() +
    theme(
        legend.position = "bottom"
    ) +
    scale_color_manual(
        values = get_gambl_colours("mutation"),
        name = "Variant Classification",
        labels = c(
            "Frame Shift In/Del",
            "Missense Mutation",
            "Nonsense Mutation"
        )
    ) +
    ggtitle("CSF2RB") +
    scale_size(guide = "none") +
    scale_fill_manual(values = c("#5c7ace", "#3baac1")) +
    theme(plot.title = element_text(face = c("bold.italic")))


ggsave("Figure_4B_CSF2RB_lollipop.pdf", width = 10, height = 4)
