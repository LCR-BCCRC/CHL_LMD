## R code to create figure 2.A in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
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
################################################
############## define custom colors ###############
################################################
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
lmd_maf <- read_tsv("../../data/Fig2/Input_Figure_2A.maf")
lmd_clin <- read_tsv("../../data/Fig2/Input_Figure_2A_clinical.tsv")

################################################
############## Create and save Oncoplot ###############
################################################

# Italicize gene names for prettyOncoplot
ht_opt(heatmap_row_names_gp = gpar(fontface = "italic"))

lmd_maf <- lmd_maf %>%
    mutate(Variant_Classification = str_replace(Variant_Classification, "Start_Codon_Del", "Translation_Start_Site")) %>%
    filter(!is.na(Variant_Classification))
lmd_clin <- lmd_clin %>%
    mutate(
        sample_id = Tumor_Sample_Barcode,
        Sex = SEX,
        EBV = case_when(
            EBERPos == 1 ~ "POS",
            EBERPos == 0 ~ "NEG"
        ),
        Mediastinal = case_when(
            Mediastinalmass == 1 ~ "POS",
            Mediastinalmass == 0 ~ "NEG"
        ),
        Subtype = DIAG_main_subsets
    )


pdf("Fig2A_oncoplot.pdf", height = 10, width = 12)
prettyOncoplot(
    these_samples_metadata = lmd_clin,
    maf = lmd_maf,
    minMutationPercent = 5,
    metadataColumns = c(
        "EBV",
        "Mediastinal",
        "Sex",
        "Subtype"
    ),
    removeNonMutated = TRUE,
    include_noncoding = list(
        "BCL7A" = c("5'UTR", "3'UTR"),
        "SOCS1" = c("5'UTR", "3'UTR")
    ),
    custom_colours = list(
        EBV = colors_EBV,
        Mediastinal = colors_mediastinal,
        Subtype = colors_DIAG,
        Sex = get_gambl_colours("sex")[c("M", "F")]
    ),
    metadataBarHeight = 4,
    metadataBarFontsize = 10,
    fontSizeGene = 10,
    legend_row = 1
)
dev.off()

# Changes made manually in Illustrator:
# - Put NA boxes in legend in black outline
# - Make combined "Frame Shift In/Del" and "In Frame In/Del" categories in legend
