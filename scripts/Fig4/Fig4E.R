## R code to create figure 4E in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
############## input data ################ 
input_Figure_4E <- read_tsv("../../data/Fig4/Fig4E_L428_CSF2RB_E788vsWT_GO_DEseq.tsv")
#############
max_dir_p = 23.82
max_x = 104
pathways_IL2 <- c("FUNG_IL2_SIGNALING_1","MARZEC_IL2_SIGNALING_UP","HALLMARK_IL2_STAT5_SIGNALING","GSE26290_WT_VS_PDK1_KO_ANTI_CD3_AND_IL2_STIM_CD8_TCELL_UP,"GSE46606_IRF4_KO_VS_WT_CD40L_IL2_IL5_1DAY_STIMULATED_BCELL_UP","GSE46606_IRF4_KO_VS_WT_CD40L_IL2_IL5_1DAY_STIMULATED_BCELL_UP)
#############
input_Figure_4E %>% 
  ggplot(aes(Enrichment,-1*log10(pvalue))) +
  geom_point(shape =21, 
             aes(size = ID %in% pathways_IL2$ID & Enrichment > 0,
                 alpha = ID %in% pathways_IL2$ID  & Enrichment > 0,
                 fill =ID %in% pathways_IL2$ID  & Enrichment > 0)) +
  coord_cartesian(xlim = c(-50,max_x)) +
  ylab("-1*log10(pvalue)") +
  geom_vline(xintercept = 0, color = "steelblue3",linetype = "dashed") +
  theme_cowplot() +
  scale_size_manual(values = c(2.5,3.2)) +
  scale_fill_manual(values = c("grey","red")) +
  scale_alpha_manual(values = c(0.02,0.8)) +
  theme(legend.position = "none") +
  ggtitle("L428")
############## save plot ################ 
ggsave("Figure4E.pdf",height = 8, width = 6)
