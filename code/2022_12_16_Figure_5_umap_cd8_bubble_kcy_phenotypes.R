# Figure 5, UMAP cluster by Color, and Bubble Map by Phenotype
require(ggplot2)
require(dplyr)
require(readr)
require(tibble)
# INPUTS:
umap   <- readr::read_tsv('data/cytof_tetramer_umap/all_cd8_umap_coordinates_and_meta_phenographs.tsv')
bubble <- readr::read_tsv('data/cytof_tetramer_umap/kcy_tet_pos_cell_by_meta_phenograph_cluster_by_pubid.csv')

### OUPUTS:
# pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_8x2.5.pdf', width = 8, height = 2.5)
# pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_5x4p5_no_legend.pdf', width = 5, height = 4.5)
# pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_7x4p5_w_legend.pdf', width = 7, height = 4.5)
# pdf('figures/2022_12_16_Figure_5_kcy_tet_bubble_plot_by_cluster.pdf', width = 3.5 , height = 4.75)


# Clusters were named by immunologists based on median marker levels,
# see data/cytof_tetramer_umap/heatmap.csv

cluster_names = c( '1'  = '1: EM',
                   '2'  = '2: TEMRA',
                   '3'  = '3: CM (endothelium homing)',
                   '4'  = '4: CM (skin homing)',
                   '5'  = '5: activated',
                   '6'  = '6: EM (resident memory)',
                   '7'  = '7: CM',
                   '8'  = '8: naïve',
                   '9'  = '9: TEMRA',
                   '10' = '10: TEMRA',
                   '11' = '11: TEMRA',
                   '12' = '12: EM',
                   '13' = '13: naïve',
                   '14' = '14: SCM',
                   '15' = '15: MAIT')

# the coarse level designation of of each cluster
cluster_coarse = c( '1'  = 'EM',
                    '2'  = 'TEMRA',
                    '3'  = 'CM',
                    '4'  = 'CM',
                    '5'  = 'activated',
                    '6'  = 'EM',
                    '7'  = 'CM',
                    '8'  = 'naïve',
                    '9'  = 'TEMRA',
                    '10' = 'TEMRA',
                    '11' = 'TEMRA',
                    '12' = 'EM',
                    '13' = 'naïve',
                    '14' = 'SCM',
                    '15' = 'MAIT')
# add cluster labels
umap = umap %>% mutate(cluster_str = as.character(MetaPhenograph)) %>% 
  mutate(cluster_name = cluster_names[cluster_str]) %>% 
  mutate(cluster_coarse = cluster_coarse[cluster_str])

fine_colors =     
  c(
    '8: naïve' =                       "#F39C12",
    '13: naïve' =                      "#FAD7A0",
    '14: SCM' =                        "#a65628",   
    '3: CM (endothelium homing)' =     "#0E6251",
    '4: CM (skin homing)' =            "#17A589",
    '7: CM' =                          "#76D7C4",
    '15: MAIT' =                       "#4daf4a",
    '5: activated' =                   "#e41a1c",
    '1: EM' =                          "#4A235A",
    '6: EM (resident memory)' =        "#7D3C98",
    '12: EM' =                         "#BB8FCE",
    '2: TEMRA' =                      '#85C1E9', 
    '9: TEMRA' =                      '#1B4F72',                   
    '10: TEMRA' =                     '#2E86C1',                    
    '11: TEMRA' =                     '#21618C')  
#https://htmlcolorcodes.com                                        
require(dplyr)
require(ggplot2)
umap
umap_copy = umap %>% select(UMAP1, UMAP2)

gg1 = ggplot(umap, aes(UMAP1, UMAP2))+ 
  geom_density_2d(data = umap_copy, color = "gray", size = .1)+
  geom_point(aes(col = cluster_name), 
             size = .001, alpha= .1) + 
  theme_bw() +
  #  facet_wrap(~MetaPhenograph) + 
  theme(strip.background = element_blank()) + 
  theme(axis.text = element_blank()) + 
  #theme(legend.position = "none") + 
  scale_color_manual("",values = fine_colors ) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

gg1  + guides(colour = guide_legend(override.aes = list(size=3, alpha = .75) ) )

pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_8x2.5.pdf', width = 8, height = 2.5)
gg1  + guides(colour = guide_legend(override.aes = list(size=3, alpha = .75))) + 
  facet_wrap(~forcats::fct_relevel(cluster_name, names(fine_colors)), 
             ncol = 8) + theme(strip.text = element_text(size = 5)) + 
  theme(legend.position = "none")
dev.off()

pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_5x4p5_no_legend.pdf', width = 5, height = 4.5)
gg1  + guides(colour = guide_legend(override.aes = list(size=3, alpha = .75))) + 
  theme(legend.position = "none")
dev.off()

pdf('figures/2022_12_16_Figure_5_cytof_umap_color_by_cluster_phenotype_7x4p5_w_legend.pdf', width = 7, height = 4.5)
gg1  + guides(colour = guide_legend(override.aes = list(size=3, alpha = .75))) + 
  theme(legend.position = "right")
dev.off()

# Square Bubble Plots
order_of_cluster_matching_above = purrr::map_chr(stringr::str_split(names(fine_colors), pattern = ":"), ~.x[1])
colors_from_above = fine_colors
names(colors_from_above) = order_of_cluster_matching_above
colors_from_above
library(ggplot2)
library(dplyr)
#df <- read.csv("bubble_plot.csv")
df1 <- bubble %>% 
  mutate(., Proportion = Cluster_freq * Tet_freq) %>% 
  left_join(tibble::tribble(~Sample,~PTID,~PTID_CODE,~Dose,
    "CO001-535869TP2"  ,"CO001-535869","P1", "0",
    "CO001-535869TP105","CO001-535869","P1", "1",
    "CO001-535869TP107","CO001-535869","P1", "3",
    "CO001-535869TP109","CO001-535869","P1", "2",
    "CO001-457348TP2"  ,"CO001-457348","P2","0",
    "CO001-457348TP104","CO001-457348","P2","1",
    "CO001-457348TP106","CO001-457348","P2","2",
    "CO001-457348TP108","CO001-457348","P2","3",
    "CO001-717584TP1"  ,"CO001-717584","P3","0",
    "CO001-717584TP103","CO001-717584","P3","1",
    "CO001-940170TP1"  ,"CO001-940170","P4","0",
    "CO001-940170TP106","CO001-940170","P4", "2",
    "CO001-282806TP1"  , "CO001-282806","XN1","0",
    "CO001-282806TP103", "CO001-282806","XN1","2",
    "CO001-282806TP104", "CO001-282806","XN1","3"))

ggdf = df1 %>% filter(Tet_type == "KCY") %>% 
  mutate(Cluster_freq = as.numeric(Cluster_freq)) %>%
  mutate(cluster_str = as.character(Cluster))

top_cluster_freq = ggdf %>% filter(Cluster_freq > .25)

pdf('figures/2022_12_16_Figure_5_kcy_tet_bubble_plot_by_cluster.pdf', width = 3.5 , height = 4.75)
ggplot(ggdf %>% filter(Cluster_freq > 0),
       aes(y= forcats::fct_rev(Dose),
           x = factor(cluster_str, levels = order_of_cluster_matching_above), 
           size  = 100*Cluster_freq,
           fill = cluster_str, col = cluster_str)) + 
  #size = log10(Proportion))) +
  geom_point(shape = 22, alpha = .7) + 
  facet_wrap(~Tet_type) + 
  theme_classic() +
  scale_fill_manual(values = colors_from_above) +
  scale_color_manual(values = colors_from_above) +
  theme(axis.text.x = element_text(angle = 0, size = 8)) + 
  facet_grid(PTID_CODE~"",scale = "free_x") + 
  theme(strip.background = element_blank()) +
  theme(strip.text.y = element_text(angle =0)) + 
  theme(legend.position = "top") + 
  theme(axis.title = element_blank()) + 
  theme(axis.text = element_text(size = 7)) +
  scale_size_continuous("", breaks = 100*c(.001, 0.01, .1, .5, 1), labels =c(0.01, 0.1, 10, 50, 100) )+
  guides(fill="none", color = "none") + 
  geom_text(data = top_cluster_freq, col ="white", aes(label = 100*round(Cluster_freq,2 )), size = 2)
dev.off()






# EXTRA
# REFERENCE COLORS 
# categorical colors
# catcolors = c("#e41a1c", # activated Red
#               "#66c2a5", # CM  teal
#               "#4daf4a", # EM - grean
#               "#984ea3", # MAIT = purple
#               "#ff7f00", # NAIVE  orang
#               "#a65628", # SCM - brown
#               "#377eb8", # TEMRA - BLUE
#               "#999999",
#               "#ffff33","#f781bf", "#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3")
#
#'57` = "CO001-535869", p1
#'80` = "CO001-457348", p2
#'44` = "CO001-717584", p3
#'37` = "CO001-940170", p4
#'04` = "CO001-282806") XN1
# coarse_colors = c('activated' = "#e41a1c", # activated RED
#                   'CM' = "#66c2a5",        # CM  teal
#                   'EM' = "#984ea3" ,        # EM -purple 
#                   'MAIT' = "#4daf4a",      # MAIT = green
#                   'naïve' = "#ff7f00",     # NAIVE  orange
#                   'SCM' = "#a65628",       # SCM - brown
#                   'TEMRA' ="#377eb8" )     # TEMRA - BLUE