# 2022_12_16_Figure_5_epitope_mapping.R

# Originally: 2022_02_01_peptides_cohen.R
# Peptide Mapping and ELISpot
# The purpose of this script is to show ELISpot experiments that led to
# confirmation of A*03 KCY immunodominatn epitope

# DEPENDENCIES: 
require(dplyr)
require(ggplot2)
require(readr)
require(ggpubr)

# INPUTS: 
#<f0> file0 contains subpool and minipool elispot data
#<f1> file1 contains individual peptide results
#<f2> contains meta-data on peptides
f0 = 'data/epitope_mapping/subpool_minipool_testing.tsv'
f1 = 'data/epitope_mapping/sars_cov2_individual_peptides.tsv'
f2 = 'data/epitope_mapping/fh_sars_cov2_spike_peptid_list.tsv'

# OUTPUTS:
figure_output_1 = 'figures/2022_12_16_Figure_5_elispot_mini_pools_no_legend.pdf'
figure_output_2 = 'figures/2022_12_16_Figure_5_elispot_peptides_in_pool_1.3.pdf'
figure_output_3 = 'figures/2022_12_16_Figure_5_elispot_peptides_in_pool_1.4.pdf'
figure_output_4 = 'figures/2022_12_16_Figure_5_elispot_mini_pools_legend.pdf'

# PREPROCESSING:
# <df_pep> - information on each peptide with clean the names
df_pep = readr::read_tsv(f2) %>% 
  select(pnum = `peptide #`,
         seq = `Sequence (COVID Spike)`, 
         pool = pool,
         sub_pool = `sub-pool`,
         mini_pool,
         number_in_subpool = `# in subpool`, 
         AA_start ,AA_end,
         S_start_end = S_start_end) 
df_pep
# <df_range> - information on the portion of S spanned by each mini pool
#   -- min_AA, starting amino acid position of the pool 
#   -- S_pool_range, e.g. S 1:67 the span of the pool 
#   We will use these for plotting, so we can be more descriptive about each pool
df_range = df_pep %>% group_by(mini_pool) %>%
  summarise(min_AA =min(AA_start),  S_pool_range =  paste0("S ", min(AA_start),":",max(AA_end))) %>% 
  mutate(mini_pool = as.character(mini_pool))
df_range
# <df0> -- dataframe with mini-pool results 
readr::read_tsv(f0) %>% pull(PUBID) %>% unique()

df0 = readr::read_tsv(f0) %>% 
  filter(PUBID %in% c("CO001-717584",
                      "CO001-535869",
                      "CO001-457348",
                      "CO001-940170")) %>%
  select(PUBID, PEPTIDESET, SFC_PER_MILLION,RESULT,REP1,REP2,	REP3) %>%
  #select(PTID, PEPTIDESET, SFC_PER_MILLION,RESULT,REP1,REP2,	REP3) %>% 
  filter(stringr::str_detect(PEPTIDESET, pattern = "COV2")) %>% 
  mutate(mini_pool = stringr::str_extract(PEPTIDESET, pattern = '[0-9]\\.[0-9]\\.?[0-9]?')) %>% 
  left_join(df_range, by = "mini_pool")

# <df_ready> Join data with peptide information
df_ready = readr::read_tsv(f1) %>% 
  select(PUBID, 
         #PTID,
         PEPTIDESET,
         AVG_BKGDSUBTRACTED,
         SFC_PER_MILLION,
         RESULT,
         REP1,REP2,REP3) %>% 
  mutate(pnum = stringr::str_extract(string = PEPTIDESET, pattern = 'CONS ([0-9]{1,3})') %>% 
           stringr::str_replace(., pattern = "CONS ",replacement="") %>%
           as.numeric()) %>% 
  left_join(df_pep, by = 'pnum') %>% 
  filter(!is.na(sub_pool)) %>% 
  mutate(RESULT = ifelse(RESULT != "Pos", "Neg", "Pos"))

# For the purposes of labeling we, subset the top peptides 
# from pool 1.3 ("S1 345:411")
df_focus_13 = df_ready %>% filter(SFC_PER_MILLION > 100) %>% 
  filter(sub_pool == "1.3") %>%
  mutate(sub_pool = "S 345:411") %>% 
  arrange(desc(SFC_PER_MILLION)) %>% 
  group_by(seq) %>% 
  slice(1)

# For the purposes of labeling we, subset the top peptides 
# from pool 1.4 ("S1 573:643")
df_focus_14 = df_ready %>% filter(SFC_PER_MILLION >100) %>% 
  filter(sub_pool == "1.4") %>%
  mutate(sub_pool = "S 573:643") %>% 
  arrange(desc(SFC_PER_MILLION)) %>% 
  group_by(seq) %>% 
  slice(1)

# GRAPHICS:
## <standard_theme>  
standard_theme = theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 7,  hjust = 1, vjust = .5 ),
        axis.text.y = element_text( size = 7 ), 
        axis.title = element_text(size = 8)) + 
  theme(legend.position = "top") + 
  theme(strip.text = element_text(face = "bold", size = 8), 
        strip.background = element_rect(colour="white", fill="#00000010")) 

## gg_peptide_13 (is the peptide level results from pool 1.3 ("S1 345:411")
gg_peptide_13 = ggplot(df_ready %>% filter(sub_pool == "1.3") %>% mutate(sub_pool = "S 345:411"), 
                       aes(x = forcats::fct_reorder(S_start_end, pnum), y =SFC_PER_MILLION)) + 
  geom_line(aes(group = factor(PUBID)), col = "gray", size = .1) + 
  geom_point(aes(shape = factor(PUBID), col= RESULT), size = .5) +
  facet_wrap(~sub_pool, scale = "free", ncol = 1) + 
  geom_text(data = df_focus_13, aes(label = seq),  size = 2, nudge_x= 0, nudge_y = 250) + 
  coord_cartesian(ylim =c(0,8000)) + 
  xlab("15-mer") +
  ylab(expression(SFC / 10^6~ PBMC)) +
  scale_shape_discrete("") + 
  scale_color_manual("",values = c("#00000020","black")) + 
  standard_theme

## gg_peptide_14 (is the peptide level results from pool 1.4 ("S1 573:643")
gg_peptide_14 = ggplot(df_ready %>% filter(sub_pool == "1.4") %>% mutate(sub_pool = "S 573:643"), 
                       aes(x = forcats::fct_reorder(S_start_end, pnum), y =SFC_PER_MILLION)) + 
  geom_line(aes(group = factor(PUBID)), col = "gray", size = .1) + 
  geom_point(aes(shape = factor(PUBID), col= RESULT), size = .5) +
  geom_text(data = df_focus_14, aes(label = seq),  size = 2, nudge_x= 0, nudge_y = 250) + 
  facet_wrap(~sub_pool, scale = "free", ncol = 1) + 
  coord_cartesian(ylim =c(0,8000)) + 
  xlab("15-mer") +
  ylab(expression(SFC / 10^6~ PBMC)) +
  scale_shape_discrete("") + 
  scale_color_manual("",values = c("#00000020","black")) + 
  standard_theme

# Mini pool results
df_line_S1 <- data.frame(x1 = 'S 1:67', x2 = 'S 573:643', y1 = 7700,y2 =  7700)
df_line_S2 <- data.frame(x1 = 'S 633:703', x2 = 'S 1133:1291', y1 =7700, y2 =  7700)
df_text_S1_S2 <- data.frame(x = c('S 345:411', 'S 889:947'), y = c(8000, 8000), label = c("S1","S2"))
df_highlight = df0 %>% filter(SFC_PER_MILLION> 2000)

gg_mini_pool = ggplot(df0 %>% 
                        mutate(header = 'S-CoV-2 Spike'), 
                      aes(x =forcats::fct_reorder(S_pool_range, min_AA), 
                          y = SFC_PER_MILLION, 
                          pch = factor(PUBID))) +
  geom_line(aes(group = factor(PUBID)), col = "gray", size = .1) + 
  geom_point(aes(col= RESULT), size = .5) + 
  facet_wrap(~header, scale = "free", ncol = 1) + 
  ylab(expression(SFC / 10^6~ PBMC)) +
  xlab("Mini-Pool") +
  coord_cartesian(ylim =c(0,8000)) + 
  scale_shape_discrete("") + 
  scale_color_manual("",values = c("#00000020","black")) + 
  standard_theme + 
  geom_segment(data = df_line_S1, 
               aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2,
                   group = NULL,
                   pch = NULL,
                   col = NULL),
               col = "gray") + 
  geom_segment(data = df_line_S2, 
               aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2,
                   group = NULL,
                   pch = NULL,
                   col = NULL),
               col = "gray") + 
  geom_text(data = df_text_S1_S2, aes(x = x, y =y, label = label,group = NULL, pch = NULL, col = NULL), col = "gray")+
  geom_text(data=df_highlight, aes(label =S_pool_range ), size = 2, nudge_x =2)

# OUTPUT FIGURES
# NAMES OF FIGURE'S SPECIFIED AT THE TOP OF SCRIPT, SHOWN HERE FOR REFERENCE
#figure_output_1 = 'figures/2022_12_16_Figure_5_elispot_mini_pools_no_legend.pdf'
#figure_output_2 = 'figures/2022_12_16_Figure_5_elispot_peptides_in_pool_1.3.pdf'
#figure_output_3 = 'figures/2022_12_16_Figure_5_elispot_peptides_in_pool_1.4.pdf'
#figure_output_4 = 'figures/2022_12_16_Figure_5_elispot_mini_pools_legend.pdf'

pdf(figure_output_1, width = 3, height=3)
gg_mini_pool + theme(legend.position = "none")
dev.off()

pdf(figure_output_2, width = 3, height=3)
gg_peptide_13 + theme(legend.position = "none")
dev.off()

pdf(figure_output_3, width = 3, height=3)
gg_peptide_14 + theme(legend.position = "none")
dev.off()

pdf(figure_output_4, width = 6, height=6)
gg_mini_pool 
dev.off()
