
###########################################  load packages ###########################################
library(data.table)
library(tidyverse)
library(pheatmap)
library(readxl)
library(ggpubr)
library(hues)
library(stringr)
library(spatstat)
library(Rphenograph)
library(openxlsx)
library(iTALK)
library(ggbeeswarm)
library(ggh4x)
library(RColorBrewer)
###########################################  Figure 3A ###########################################
# CST
p.data <- read.csv(paste0('data/Fig3/CST_patient.csv'))
location_X <- c(1000, 1500)
location_Y <- c(0, 500)

## visualizaiton of the whole image
ggplot(p.data, aes(x = Location_Center_X, y = Location_Center_Y)) +
  geom_point(data = p.data %>% filter(initial_cell_type == 'HRS'), color = 'red', size = 2) +
  geom_point(data = p.data %>% filter(FOXP3_CD4 == 1), color = '#33a02c', size = 2) +
  geom_point(data = p.data %>% filter(LAG3_CD4 == 1), color = '#FFD100', size = 2) +
  geom_point(data = p.data %>% filter(GranzymeB_CD8 == 1), color = '#710280', size = 2) +
  geom_point(data = p.data %>% filter(CD163_Mac_Myeloid == 1), color = '#fb9a99', size = 2) +
  geom_point(data = p.data %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black')) +
  annotate("rect", xmin = location_X[1], xmax = location_X[2], 
           ymin = location_Y[1], ymax = location_Y[2],
           alpha = 0, color = "white", lwd = 3) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -40, ymax = -40,
           alpha = 0, color = "white", lwd = 3) +
  annotate("text", x = 50, y = -20, fontface = "bold", size = 7, color = 'white',
           label = '100 μm')


## visualization of the enriched region
p.data.subset <- p.data %>%
  filter(Location_Center_X > location_X[1],
         Location_Center_X < location_X[2],
         Location_Center_Y > location_Y[1],
         Location_Center_Y < location_Y[2])
ggplot(p.data.subset) +
  geom_point(data = p.data.subset%>% filter(initial_cell_type == 'HRS'), aes(x = Location_Center_X, y = Location_Center_Y), color = 'red', size = 2) +
  geom_point(data = p.data.subset %>% filter(FOXP3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#33a02c', size = 2) +
  geom_point(data = p.data.subset %>% filter(LAG3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#CAB651', size = 2) +
  geom_point(data = p.data.subset %>% filter(GranzymeB_CD8 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#6a3d9a', size = 2) +
  geom_point(data = p.data.subset %>% filter(CD163_Mac_Myeloid == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#fb9a99', size = 2) +
  geom_point(data = p.data.subset %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))+
  annotate("rect", xmin = location_X[1] , xmax = location_X[1] + 50, ymin = location_Y[1] - 20, ymax = location_Y[1] -20,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 25, y = location_Y[1] - 8, fontface = "bold", size = 8, color = 'white',
           label = '50 μm')


# CN913
p.data <- read.csv(paste0('data/Fig3/CN913_patient.csv'))
location_X <- c(200, 700)
location_Y <- c(500, 1000)

## visualizaiton of the whole image
ggplot(p.data, aes(x = Location_Center_X, y = Location_Center_Y)) +
  geom_point(data = p.data %>% filter(initial_cell_type == 'HRS'), color = 'red', size = 2) +
  geom_point(data = p.data %>% filter(FOXP3_CD4 == 1), color = '#33a02c', size = 2) +
  geom_point(data = p.data %>% filter(LAG3_CD4 == 1), color = '#FFD100', size = 2) +
  geom_point(data = p.data %>% filter(GranzymeB_CD8 == 1), color = '#710280', size = 2) +
  geom_point(data = p.data %>% filter(CD163_Mac_Myeloid == 1), color = '#fb9a99', size = 2) +
  geom_point(data = p.data %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black')) +
  annotate("rect", xmin = location_X[1], xmax = location_X[2], 
           ymin = location_Y[1], ymax = location_Y[2],
           alpha = 0, color = "white", lwd = 3) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -40, ymax = -40,
           alpha = 0, color = "white", lwd = 3) +
  annotate("text", x = 50, y = -20, fontface = "bold", size = 7, color = 'white',
           label = '100 μm')


## visualization of the enriched region
p.data.subset <- p.data %>%
  filter(Location_Center_X > location_X[1],
         Location_Center_X < location_X[2],
         Location_Center_Y > location_Y[1],
         Location_Center_Y < location_Y[2])
ggplot(p.data.subset) +
  geom_point(data = p.data.subset%>% filter(initial_cell_type == 'HRS'), aes(x = Location_Center_X, y = Location_Center_Y), color = 'red', size = 2) +
  geom_point(data = p.data.subset %>% filter(FOXP3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#33a02c', size = 2) +
  geom_point(data = p.data.subset %>% filter(LAG3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#CAB651', size = 2) +
  geom_point(data = p.data.subset %>% filter(GranzymeB_CD8 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#6a3d9a', size = 2) +
  geom_point(data = p.data.subset %>% filter(CD163_Mac_Myeloid == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#fb9a99', size = 2) +
  geom_point(data = p.data.subset %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))+
  annotate("rect", xmin = location_X[1] , xmax = location_X[1] + 50, ymin = location_Y[1] - 20, ymax = location_Y[1] -20,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 25, y = location_Y[1] - 8, fontface = "bold", size = 8, color = 'white',
           label = '50 μm')


# STB
p.data <- read.csv(paste0('data/Fig3/STB_patient.csv'))
location_X <- c(400, 900)
location_Y <- c(700, 1200)

## visualizaiton of the whole image
ggplot(p.data, aes(x = Location_Center_X, y = Location_Center_Y)) +
  geom_point(data = p.data %>% filter(initial_cell_type == 'HRS'), color = 'red', size = 2) +
  geom_point(data = p.data %>% filter(FOXP3_CD4 == 1), color = '#33a02c', size = 2) +
  geom_point(data = p.data %>% filter(LAG3_CD4 == 1), color = '#FFD100', size = 2) +
  geom_point(data = p.data %>% filter(GranzymeB_CD8 == 1), color = '#710280', size = 2) +
  geom_point(data = p.data %>% filter(CD163_Mac_Myeloid == 1), color = '#fb9a99', size = 2) +
  geom_point(data = p.data %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black')) +
  annotate("rect", xmin = location_X[1], xmax = location_X[2], 
           ymin = location_Y[1], ymax = location_Y[2],
           alpha = 0, color = "white", lwd = 3) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -40, ymax = -40,
           alpha = 0, color = "white", lwd = 3) +
  annotate("text", x = 50, y = -20, fontface = "bold", size = 7, color = 'white',
           label = '100 μm')


## visualization of the enriched region
p.data.subset <- p.data %>%
  filter(Location_Center_X > location_X[1],
         Location_Center_X < location_X[2],
         Location_Center_Y > location_Y[1],
         Location_Center_Y < location_Y[2])
ggplot(p.data.subset) +
  geom_point(data = p.data.subset%>% filter(initial_cell_type == 'HRS'), aes(x = Location_Center_X, y = Location_Center_Y), color = 'red', size = 2) +
  geom_point(data = p.data.subset %>% filter(FOXP3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#33a02c', size = 2) +
  geom_point(data = p.data.subset %>% filter(LAG3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#CAB651', size = 2) +
  geom_point(data = p.data.subset %>% filter(GranzymeB_CD8 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#6a3d9a', size = 2) +
  geom_point(data = p.data.subset %>% filter(CD163_Mac_Myeloid == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#fb9a99', size = 2) +
  geom_point(data = p.data.subset %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))+
  annotate("rect", xmin = location_X[1] , xmax = location_X[1] + 50, ymin = location_Y[1] - 20, ymax = location_Y[1] -20,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 25, y = location_Y[1] - 8, fontface = "bold", size = 8, color = 'white',
           label = '50 μm')


# CN2P
p.data <- read.csv(paste0('data/Fig3/CN2P_patient.csv'))

location_X <- c(200, 700)
location_Y <- c(1000, 1500)

## visualizaiton of the whole image
ggplot(p.data, aes(x = Location_Center_X, y = Location_Center_Y)) +
  geom_point(data = p.data %>% filter(initial_cell_type == 'HRS'), color = 'red', size = 2) +
  geom_point(data = p.data %>% filter(FOXP3_CD4 == 1), color = '#33a02c', size = 2) +
  geom_point(data = p.data %>% filter(LAG3_CD4 == 1), color = '#FFD100', size = 2) +
  geom_point(data = p.data %>% filter(GranzymeB_CD8 == 1), color = '#710280', size = 2) +
  geom_point(data = p.data %>% filter(CD163_Mac_Myeloid == 1), color = '#fb9a99', size = 2) +
  geom_point(data = p.data %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black')) +
  annotate("rect", xmin = location_X[1], xmax = location_X[2], 
           ymin = location_Y[1], ymax = location_Y[2],
           alpha = 0, color = "white", lwd = 3) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -40, ymax = -40,
           alpha = 0, color = "white", lwd = 3) +
  annotate("text", x = 50, y = -20, fontface = "bold", size = 7, color = 'white',
           label = '100 μm')

## visualization of the enriched region
p.data.subset <- p.data %>%
  filter(Location_Center_X > location_X[1],
         Location_Center_X < location_X[2],
         Location_Center_Y > location_Y[1],
         Location_Center_Y < location_Y[2])

ggplot(p.data.subset) +
  geom_point(data = p.data.subset%>% filter(initial_cell_type == 'HRS'), aes(x = Location_Center_X, y = Location_Center_Y), color = 'red', size = 2) +
  geom_point(data = p.data.subset %>% filter(FOXP3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#33a02c', size = 2) +
  geom_point(data = p.data.subset %>% filter(LAG3_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#CAB651', size = 2) +
  geom_point(data = p.data.subset %>% filter(GranzymeB_CD8 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#6a3d9a', size = 2) +
  geom_point(data = p.data.subset %>% filter(CD163_Mac_Myeloid == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#fb9a99', size = 2) +
  geom_point(data = p.data.subset %>% filter(PD1_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#095f80', size = 2) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))+
  annotate("rect", xmin = location_X[1] , xmax = location_X[1] + 50, ymin = location_Y[1] - 20, ymax = location_Y[1] -20,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 25, y = location_Y[1] - 8, fontface = "bold", size = 8, color = 'white',
           label = '50 μm')


###########################################  Figure 3B ###########################################
# load data
expr.dt <- read.csv(paste0('data/Fig3/expr_dt.csv'))

# order cell types and NMF clusters
expr.dt$cell_type <- factor(expr.dt$cell_type, levels = c('PDL1+ HRS', 'TARC+ HRS', 'PDL1+ Mac', 'CD68+ Mac', 'CD163+ Mac', 'CTL', 'PD1+ CD4', 'TFH','LAG3+ CD4','FOXP3+ CD4', 'CD27+ B'))
expr.dt$Nlee4 <- factor(expr.dt$Nlee4, levels = c('CST', 'CN913', 'STB', 'CN2P'))

# visualization
ggplot(expr.dt) +
  geom_point(aes(x = Nlee4, y = cell_type, fill = scaled, size = log_p_adj), shape = 21, color = 'black' ) +
  labs(x = '', y ='', fill = 'Spatial score (scaled)', size = '-Log10(p-value)')+
  theme_classic() +
  facet_grid2(.~Nlee4, scales = "free_x", space = "free_x", 
              remove_labels = "all",
              strip = strip_themed(
                background_x = elem_list_rect(
                  fill = c('#A6FFB6', '#33A1F6', '#D16303', '#9370DB')),
                text_x = element_text(size = 80))) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 80),
        axis.title = element_text(size = 80),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 50),
        panel.spacing.x  = unit(0, "mm"),
        legend.text = element_text(size = 50))+
  scale_fill_gradientn(colours = c('#00007F','blue','#007FFF','lightblue','white',
                                   'orange', '#FF7F00', 'red','#7F0000')) +
  scale_size_area(max_size = 60) +
  ggpubr::border()



###########################################  Figure 3C ###########################################
# load data
cellchat <- read.csv(file = paste0('data/Fig3/cellchat.csv'))

# order NMF clusters
cellchat$cluster <- factor(cellchat$cluster,
                           levels = rev(c('CN2P',
                                          'STB',
                                          'CN913',
                                          'CST')))

# order cell interaction pairs
cellchat$cell_pair <- factor(cellchat$cell_pair, levels = rev(c('TARC+ HRS > CCR4+ CD4',
                                                                "TARC+ HRS > FOXP3+ CD4",
                                                                "CXCL13+ CD4 > CXCR5+ HRS",
                                                                "PDL1+ HRS > PD1+ CD4",
                                                                "PDL1+ HRS > TFH")))

# order ligand-receptor pairs
cellchat$interaction_name <- factor(cellchat$interaction_name, levels = c('TARC - CCR4', 'CXCL13 - CXCR5', 'PDL1 - PD1'))

strip <- strip_themed(background_x = elem_list_rect(fill = c('#A6FFB6', '#33A1F6', '#D16303', '#9370DB')))

ggplot(cellchat, aes( x= interaction_name, y = cell_pair, fill = interaction_strength)) +
  geom_tile(size = 8) +
  facet_grid2(~cluster,strip = strip)+
  theme(axis.text.x = element_text(size = 60, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 60),
        axis.title = element_text(size = 60),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 40),
        legend.key.size = unit(2, 'cm'),
        strip.text = element_text(size = 60)) +
  labs(x = '', y = '', fill = 'Interaction strength') +
  ggpubr::border() +
  scale_fill_gradientn(colors= colorRampPalette(rev(brewer.pal(n = 7, 
                                                               name = "RdYlBu")))(100)) 



###########################################  Figure 3D ###########################################
# color code
cell_col <- c('TARC+ HRS' = '#e31a1c',
              'CCR4+ CD4' = '#b3de69',
              'Galectin9+ HRS' = "#C8539B",
              'LAG3+ CD4' = '#CAB651',
              'PD1+ CD4' = '#1f78b4',
              'PDL1+ HRS' = "#78CF55",
              'CD80+ HRS' = "#C0523B",
              'CTLA4+ CD4' = "#4B2C46",
              'CXCL13+ CD4' = "#85C9B3",
              'FOXP3+ CD4'= "#33a02c")

# load data
circus.dt <- read.csv(paste0('data/Fig3/circus_dt.csv'))

# visualization
LRPlot(circus.dt,datatype='mean count',link.arr.lwd=circus.dt$mean_expr,link.arr.col = circus.dt$line_col, cell_col = cell_col,transparency = 0.8)



###########################################  Figure 3E ###########################################
# load data
CCR4.patient.spatial.score <- read.csv(paste0('data/Fig3/CCR4_patient_spatial_score.csv'))

# visualization
ggplot(CCR4.patient.spatial.score,aes(x = CCL17_status, y= mean)) +
  geom_boxplot(aes(x = CCL17_status, y= mean, color = CCL17_status), linewidth = 3)+
  geom_quasirandom(data = CCR4.patient.spatial.score,aes(x = CCL17_status, y= mean, color = CCL17_status), size = 6) +
  stat_compare_means(method = "t.test", comparisons = list(c("TARC+ HRS", "TARC- HRS")), 
                     paired = FALSE, label.x.npc = 0.5, hjust = 0.5, size = 20, tip.length = 0,label = "p.signif")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 80, colour = 'black'),
        axis.text.y = element_text(size = 50),
        axis.title = element_text(size = 80),
        legend.position = 'none') +
  labs(x ='', y = 'CCR4+ T Spatial Score')+ 
  scale_color_manual(values = c('TARC+ HRS' = '#e31a1c',
                                'TARC- HRS' = '#80b1d3')) +
  scale_fill_manual(values = c('TARC+ HRS' = '#e31a1c',
                               'TARC- HRS' = '#80b1d3')) +
  ggpubr::border()

###########################################  Figure 3F ###########################################
# load data
TARC.neg.patient <- read.csv(paste0('data/Fig3/TARC_neg_patient.csv'))
location_X <- c(900, 1650)
location_Y <- c(900, 1650)

ggplot(TARC.neg.patient) +
  geom_point(data = TARC.neg.patient %>% filter(CXCR5_HRS ==1 | CD47_HRS ==1 | CXCR4_HRS == 1 | TGFB_HRS == 1 | Galectin9_HRS ==1 | PDL1_HRS ==1 | CLEC2D_HRS ==1 | CD123_HRS ==1 | CD80_HRS ==1| Ki67_HRS == 1 & CCL17_HRS ==0), aes(x = Location_Center_X, y = Location_Center_Y), color = '#80b1d3', size = 3) +
  geom_point(data =TARC.neg.patient %>% filter(CCR4_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#b3de69', size = 3) +  
  annotate("rect", xmin = location_X[1] + 50 , xmax = location_X[1] + 80, ymin = location_Y[1] + 50, ymax = location_Y[1] +50,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 65, y = location_Y[1] + 67, fontface = "bold", size = 6, color = 'white',
           label = '30 μm') +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))

# load data
TARC.pos.patient <- read.csv(paste0('data/Fig3/TARC_pos_patient.csv'))
location_X <- c(900, 1650)
location_Y <- c(900, 1650)

ggplot(TARC.pos.patient ) +
  geom_point(data = TARC.pos.patient  %>% filter(CCL17_HRS == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#e31a1c', size = 3) +
  geom_point(data = TARC.pos.patient  %>% filter(CCR4_CD4 == 1), aes(x = Location_Center_X, y = Location_Center_Y), color = '#b3de69', size = 3) +
  annotate("rect", xmin = location_X[1] + 50 , xmax = location_X[1] + 80, ymin = location_Y[1] + 50, ymax = location_Y[1] +50,
           alpha = 0,
           color = "white",
           lwd = 3) +
  annotate("text", x = location_X[1] + 65, y = location_Y[1] + 67, fontface = "bold", size = 6, color = 'white',
           label = '30 μm') +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'black'))


