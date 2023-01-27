# 2022_12_16_Figure_2_ics_corrlation

require(dplyr)
require(readr)
require(ggplot2)
require(corrplot)
require(tidyr)

# INPUT:
request_df2 = readr::read_tsv('data/ics_correlation/ics_basline_post1_correlation_data.tsv')

# OUTPUT:
figure_output_1 = "figures/2022_12_16_Figure_2_cross_plots_priority_4_quadrant.pdf"
figure_output_2 = 'figures/2022_12_16_Figure_2_spearman_correlation_matrix.pdf'

varx = "PRE CD8+ S1"
vary = "POST CD8+ S1"
u = request_df2 %>% filter(shortname == varx ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
v = request_df2 %>% filter(shortname == vary ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
duv = left_join(u,v, by = c("pubid") ) %>%
  mutate(fake_group = 1) 
duv
ggplot(data  = duv , aes(x = pctpos.x, 
                         y = pctpos.y, 
                         shape = factor(response.mimosa.y))) + 
  geom_point() + scale_y_log10()+ scale_x_log10() + 
  theme_classic() + annotation_logticks(size = .2) + 
  scale_shape_manual("Responder", values = c(19,19)) + 
  theme(legend.position = "none") +  
  ylab(vary) + xlab(varx) + 
  #geom_point(data = filter(duv, ND_pre.x == 1), col = "black", pch =3, size =3) +
  coord_cartesian(xlim = c(0.001, .5), ylim = c(0.001, .5)) + 
  geom_abline(col = "gray", linetype = "dashed", size = .2) + 
  geom_smooth(aes(group =fake_group), method = "lm", se = F, col = "gray", size = .2) + 
  ggpubr::stat_cor(method = "spearman", aes(group = fake_group), label.x.npc =.25, label.y.npc =0.1, size = 2.5) -> gg1

gg1

varx = "PRE CD4+ S1"
vary = "POST CD8+ S1"
u = request_df2 %>% filter(shortname == varx ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
v = request_df2 %>% filter(shortname == vary ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
duv = left_join(u,v, by = c("pubid") ) %>%
  mutate(fake_group = 1) 

ggplot(data  = duv , aes(x = pctpos.x, y = pctpos.y, shape = factor(response.mimosa.y))) + 
  geom_point() + scale_y_log10()+ scale_x_log10() + 
  theme_classic() + annotation_logticks(size = .2) + 
  scale_shape_manual("Responder", values = c(19,19)) + 
  theme(legend.position = "none") + 
  ylab(vary) + xlab(varx) + 
  #geom_point(data = filter(duv, ND_pre.x == 1), col = "black", pch =3, size = 3) +
  coord_cartesian(xlim = c(0.001, .5), ylim = c(0.001, .5)) + 
  #geom_abline(col = "gray", linetype = "dashed", size = .2) + 
  geom_smooth(aes(group =fake_group), method = "lm", se = F, col = "gray", size = .2) + 
  ggpubr::stat_cor(method = "spearman", aes(group = fake_group), label.x.npc =.25, label.y.npc =0.1, size = 2.5) -> gg2

varx = "PRE CD4+ S1"
vary = "POST RBD IgG"
u = request_df2 %>% filter(shortname == varx ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
v = request_df2 %>% filter(shortname == vary ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
duv = left_join(u,v, by = c("pubid") ) %>%
  mutate(fake_group = 1) 

ggplot(data  = duv , aes(x = pctpos.x, y = pctpos.y)) + 
  geom_point( pch = 19) + scale_y_log10()+ scale_x_log10() + 
  theme_classic() + annotation_logticks(size = .2) + 
  scale_shape_manual("Responder", values = c(19,19)) + 
  theme(legend.position = "none") + 
  ylab(vary) + xlab(varx) + 
  #geom_point(data = filter(duv, ND_pre.x == 1), col = "black", pch =19, size = 1) +
  coord_cartesian(xlim = c(0.001, .5)) + 
  geom_abline(col = "gray", linetype = "dashed", size = .2) + 
  geom_smooth(aes(group =fake_group), method = "lm", se = F, col = "gray", size = .2) + 
  ggpubr::stat_cor(method = "spearman", aes(group = fake_group), label.x.npc =.25, label.y.npc =0.1, size = 2.5) -> gg3


varx = "PRE CD8+ S1"
vary = "POST RBD IgG"
u = request_df2 %>% filter(shortname == varx ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
v = request_df2 %>% filter(shortname == vary ) %>% 
  mutate(ND_pre = ifelse(pctpos < 0.001, 1, 0)) %>%
  mutate(pctpos = ifelse(pctpos < 0.001, 0.001, pctpos))
duv = left_join(u,v, by = c("pubid") ) %>%
  mutate(fake_group = 1) 

ggplot(data  = duv , aes(x = pctpos.x, y = pctpos.y)) + 
  geom_point( pch = 19) + scale_y_log10()+ scale_x_log10() + 
  theme_classic() + annotation_logticks(size = .2) + 
  scale_shape_manual("Responder", values = c(19,19)) + 
  theme(legend.position = "none") + 
  ylab(vary) + xlab(varx) + 
  #geom_point(data = filter(duv, ND_pre.x == 1), col = "black", pch =19, size = 1) +
  coord_cartesian(xlim = c(0.001, .5)) + 
  geom_abline(col = "gray", linetype = "dashed", size = .2) + 
  geom_smooth(aes(group =fake_group), method = "lm", se = F, col = "gray", size = .2) + 
  ggpubr::stat_cor(method = "spearman", aes(group = fake_group), label.x.npc =.25, label.y.npc =0.1, size = 2.5) -> gg4

#figure_output_1 = "figures/cross_plots_priority_4_quadrant.pdf",
pdf(figure_output_1, width = 5, height = 5)
gridExtra::grid.arrange(gg2, gg1, gg3, gg4)
dev.off()


request_df2 %>% group_by(pubid, visitno, shortname) %>% 
  summarise(pctpos = median(pctpos))%>%
  ungroup()%>%
  select(-visitno) %>% 
  tidyr::spread(key = shortname , value = pctpos) -> request_wide 

request_wide
#pdf("figures/Explode.pdf", width =8, height = 8)
rownames(request_wide)
request_wide %>% 
  select(-pubid) %>%
  as.matrix() %>% 
  cor(method = "spearman", use = "pairwise.complete") -> M

pdf(figure_output_2,  width =8, height = 8)
corrplot::corrplot(M, 
                     addCoef.col="black",
                     #type = "upper",
                     number.cex=.5,
                     rect.col = "white",
                     tl.col = 'black',
                     tl.srt=45,
                     cl.ratio = 0.2,  
                     diag = F,
                     tl.cex=.75) 
dev.off()
