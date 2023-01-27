# 2022_12_16_Figure_5_umap_tetramer_grid.R

# NOTES: originally heeju/20221128_Figure5.R refactored to use PUBID instead of ptids

# DEPENDENCIES
require(ggplot)
require(readr)
require(forcats)
require(dplyr)

# INPUTS
# ALL CD8+ sorted cells
dfallcd8 = readr::read_tsv('data/cytof_tetramer_umap/umap_coordinates_KCYtet_all_cd8_cells_20221128.tsv')
# ALL CD8+ KCY+ cells 
dfall    = readr::read_tsv('data/cytof_tetramer_umap/umap_coordinates_KCYtet_pos_20221128.tsv')
# Sampling timing for the 5 ptids in above plots
meta     = readr::read_tsv('data/cytof_tetramer_umap/umap_coordinates_pubid_metadata.tsv') 

# ORDER OR SAMPLES, P1, P2, P3, P4, Naive
spec_levels = c(
  "CO001-535869", 
  "CO001-457348", 
  "CO001-717584", 
  "CO001-940170",
  "CO001-282806")

dfall %>% filter(pubid == "CO001-535869") %>% group_by(visit, dose) %>% slice(1)

# OUTPUTS (A large PDF with many point 24.6 MB)
output_figure = 'figures/2022_12_16_Figure_5_umap_grid_tetramer_kcy_over_all_CD8_20221128.pdf'
output_figure_png = 'figures/2022_12_16_Figure_5_umap_grid_tetramer_kcy_over_all_CD8_20221128.png'

# SET FACTOR LEVELS
dfall$pubid = factor(dfall$pubid, levels = spec_levels)
dfallcd8pubid = factor(dfallcd8$pubid, levels = spec_levels)

pdf(output_figure , height = 8, width = 8)
ggplot(dfall %>% filter(dose != 1.9)) + 
  geom_point(data =dfallcd8, aes(x= UMAP1, y = UMAP2, alpha = NULL), alpha=.1, size = .001, col = "gray") +
  geom_point(aes(x= UMAP1, y = UMAP2, alpha = emphasis, size = emphasis), col = "red", pch = 20) +
  facet_grid(forcats::fct_relevel(pubid, spec_levels)~dose) + 
  scale_alpha_manual(values = c(.05, .5))+
  scale_size_manual(values = c(.001, .001))+
  theme_bw() +
  theme(strip.background = element_blank()) + 
  theme(axis.text = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = "None")
dev.off()

png(output_figure_png, height = 1000, width = 1000)
ggplot(dfall %>% filter(dose != 1.9)) + 
  geom_point(data =dfallcd8, aes(x= UMAP1, y = UMAP2, alpha = NULL), alpha=.1, size = .001, col = "gray") +
  geom_point(aes(x= UMAP1, y = UMAP2, alpha = emphasis, size = emphasis), col = "red", pch = 20) +
  facet_grid(forcats::fct_relevel(pubid, spec_levels)~dose) + 
  scale_alpha_manual(values = c(.05, .5))+
  scale_size_manual(values = c(.001, .001))+
  theme_bw() +
  theme(strip.background = element_blank()) + 
  theme(axis.text = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = "None")
dev.off()

meta
