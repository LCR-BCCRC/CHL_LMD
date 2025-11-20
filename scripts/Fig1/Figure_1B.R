## R code to create figure 1B in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
######################################################### 
############ load packages ############
######################################################### 
library(readr)
library(dplyr)
library(magrittr)
library(stats)
library(ggplot2)
library(cowplot)
################################################
source("get_gambl_colours.R")
source("create_onco_matrix.R")
source("prettyOncoplot.R")
################################################
############## define custom colors ###############
################################################
# Colours
colors_MHC <- c(
    "POS" = "blue4",
    "NEG" = "lightgrey"
)

colors_mediastinal <- c(
    "POS" = "firebrick2",
    "NEG" = "lightgrey",
    "NA" = "white"
)
colors_EBV <- c(
    "POS" = "chartreuse4",
    "NEG" = "lightgrey",
    "NA" = "white"
)

colors_AGE <- c(
    "<45yo" = "lightgrey",
    ">=45yo" = "cyan3"
)
colors_DIAG <- c(
    "NS" = "steelblue4",
    "MC" = "#1C871C",
    "Other" = "lightgrey"
)

colors_maf5 <- c(
    "Yes" = "firebrick2",
    "No" = "lightgrey"
)

colors_LDA <- c(
    "H1" = "red",
    "H2" = "steelblue"
)

hlgen_cols <- c(
    CST = "#A6FFB6",
    CN913 = "#33A1F6",
    STB = "#D16303",
    CN2P = "#9370DB",
    Unclassified = "#D3D3D3"
)
################################################
############## Input files ###############
################################################
#########################################################
############ Input data ############
#########################################################


# Load LMD HLGen calls and mutations
lmd_hlgen <- read_tsv("../../data/Fig1/Input_Figure_1B.tsv") %>%
    rename(cluster = HLGen) %>%
    left_join(hlgen_key) %>%
    mutate(HLGen = factor(HLGen, levels = names(hlgen_cols)[1:4]))

# Convert mutation status to matrix
hlgen_mat <- lmd_hlgen %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(matches("_MUT|_DEL|_AMP"), ~ case_when(
        str_detect(cur_column(), "DEL") ~ str_replace(str_replace(.x, "1", "Deletion"), "0", ""),
        str_detect(cur_column(), "AMP") ~ str_replace(str_replace(.x, "1", "Amplification"), "0", ""),
        str_detect(cur_column(), "MUT") ~ str_replace(str_replace(.x, "1", "Mutation"), "0", "")
    ))) %>%
    rename_with(~ str_remove(.x, "_MUT|CNA_")) %>%
    select(res_id, any_of(unname(unlist(hlgen_genes)))) %>%
    column_to_rownames("res_id") %>%
    t()

# Set column split
column_split <- lmd_hlgen %>%
    select(res_id, HLGen) %>%
    column_to_rownames("res_id")
column_split <- column_split[colnames(hlgen_mat), , drop = FALSE]

# Calculate proprtion of mutations in each cluster
props <- lmd_hlgen %>%
    rename_with(~ str_remove(.x, "_MUT|CNA_")) %>%
    select(all_of(rownames(hlgen_mat)), HLGen) %>%
    pivot_longer(-HLGen) %>%
    group_by(name) %>%
    mutate(num_mutated = sum(value)) %>%
    group_by(name, HLGen) %>%
    mutate(num_this_group = sum(value)) %>%
    select(-value) %>%
    distinct() %>%
    mutate(prop_this_group = num_this_group / num_mutated) %>%
    ungroup() %>%
    select(HLGen, name, prop_this_group) %>%
    pivot_wider(
        names_from = "HLGen",
        values_from = "prop_this_group"
    ) %>%
    column_to_rownames("name")

props <- props[, names(hlgen_cols)[1:4]]

# Create row annotation with proportions
ra <- rowAnnotation(Proportion = anno_barplot(props,
    gp = gpar(fill = hlgen_cols),
    bar_width = 1, width = unit(2, "cm")
))
row_split <- factor(unlist(lapply(c(1:4), function(x) rep(names(hlgen_genes)[x], length(hlgen_genes[[x]])))), levels = names(hlgen_cols))

# Create column annotations with HLGen classes
ha <- HeatmapAnnotation(
    HLGen = column_split$HLGen,
    col = list(HLGen = hlgen_cols),
    show_legend = FALSE
)

# EBV status, thymic mass, age, and surface MHC-I and -II

annot_bottom <- lmd_hlgen %>%
    select(
        res_id,
        EBV = EBERPos,
        Mediastinal = Mediastinalmass,
        Age = AGE_group45,
        MHCI = MHCIPOS,
        MHCII = MHCIIPOS,
        HLGen
    ) %>%
    mutate(
        across(matches("MHC|EBV|Mediastinal"), ~ case_when(
            .x == 1 ~ "POS",
            .x == 0 ~ "NEG"
        )),
    ) %>%
    mutate(Age = ifelse(Age == 1, "<45yo", ">=45yo")) %>%
    column_to_rownames("res_id")

annot_bottom <- annot_bottom[rownames(column_split), ]

ha_bottom <- HeatmapAnnotation(
    EBV = annot_bottom$EBV,
    Mediastinal = annot_bottom$Mediastinal,
    Age = annot_bottom$Age,
    MHCI = annot_bottom$MHCI,
    MHCII = annot_bottom$MHCII,
    col = list(
        EBV = colors_EBV,
        Mediastinal = colors_mediastinal,
        Age = colors_AGE,
        MHCI = colors_MHC,
        MHCII = colors_MHC
    ),
    show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    na_col = "white"
)

# Set colour scheme for matrix/heatmap
col <- c(Amplification = "firebrick", Deletion = "darkblue", Mutation = "black")

# OncoPrint
op <- oncoPrint(hlgen_mat,
    alter_fun = list(
        background = function(x, y, w, h) {
            grid.rect(x, y, w, h,
                gp = gpar(fill = "#e9e8e8", col = NA)
            )
        },
        Amplification = function(x, y, w, h) {
            grid.rect(x, y, w, h,
                gp = gpar(fill = col["Amplification"], col = NA)
            )
        },
        Deletion = function(x, y, w, h) {
            grid.rect(x, y, w, h,
                gp = gpar(fill = col["Deletion"], col = NA)
            )
        },
        Mutation = function(x, y, w, h) {
            grid.rect(x, y, w, h,
                gp = gpar(fill = col["Mutation"], col = NA)
            )
        }
    ), col = col,
    left_annotation = ra,
    bottom_annotation = ha_bottom,
    column_split = column_split$HLGen,
    row_split = row_split,
    row_order = rownames(hlgen_mat),
    cluster_row_slices = FALSE,
    top_annotation = ha,
    right_annotation = NULL,
    show_pct = FALSE
)

pdf("LMD_HLGen_oncoprint.pdf", height = 9, width = 12)
draw(op, padding = unit(c(2, 2, 2, 10), "mm"), annotation_legend_side = "bottom")
dev.off()
