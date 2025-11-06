## R code to create figure 2.B in Aoki et al, Multidimensional characterization of cellular ecosystems in Hodgkin lymphoma
############## load packages ################ 
library(tidyverse)
library(ggplot2)
library(ggsankey)
############## define colors ################ 
colors_HLGen <- c("CST" = "#A6FFB6",
                  "CN913" = "#33A1F6",
                  "STB" = "#D16303",
                  "CN2P" = "#9370DB",
                  "UNCLASS" = "#D3D3D3")

colors_custom <- c(colors_HLGen,"OLD" = "#D3D3D3","YOUNG" = "cyan3","NEG" = "#D3D3D3","POS" = "#458B00")

############## input data ################ 
input_Figure_2B <- read_tsv("../../data/input_Figure_2B_EBV_HLGen.tsv")
#############
df_sankey <- input_Figure_2B
df_sankey_long <- df_sankey %>% 
  make_long(EBV,HLGen)
df_sankey_long$node <- factor(df_sankey_long$node,levels= c("NEG","POS","CN2P","STB","CN913","CST")) ## Change order
df_sankey_long %>%
  ggplot(aes(x = x
             , next_x = next_x
             , node = node
             , next_node = next_node
             , fill = factor(node)
             , label = node)
  ) +
  geom_sankey(space = 0,
              flow.alpha = 0.5
              , node.color = "black"
              ,show.legend = FALSE,
              width = .4,
  ) +
  scale_x_discrete(position = "top") + # put axis (labels) on top of the graph
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 18,face = "bold",vjust =-2),
        panel.grid = element_blank(),
        panel.background = element_blank())  +
  scale_fill_manual(values = colors_custom_adj) +
  geom_sankey_text(space = 0,size = 5,type = "sankey") +## space = 0 to match  
  theme(axis.text.x.top = element_text(margin = margin(b = -20))) # Solution label position## space = 0 to match  
############## save plot ################ 
ggsave("Figure2B_EBV_HLGen.pdf",height = 8, width = 6)
