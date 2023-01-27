# 2022_12_16_Figure_1_ics_metaclustered.R
# FIGURE_1 Plot
# This was refactored to use PUBIDs
# Taken directly from (~/active/mcelrath_manucript/mansucript_prep/pdf_v4/parts/fig1_plot_script_10_27_2022.R)

# INPUTS:
path1 = "data/ics_metaclustered/"
x_cd8_s = readr::read_csv(paste0(path1,'ics_cd8_s.csv'))
x_cd4_s = readr::read_csv(paste0(path1,'ics_cd4_s.csv'))
x_cd8_s_subset = readr::read_csv(paste0(path1,'ics_cd8_s_N_CM_EM_TEMRA.csv'))
x_cd4_s_subset = readr::read_csv(paste0(path1,'ics_cd4_s_N_CM_EM_TEMRA.csv'))

# OUTPUTS:
figure_output_1 = 'figures/2022_12_16_Figure_1_cd4_cd8_S_INFg_IL2.pdf'
figure_output_2 = "figures/2022_12_16_Figure_1_cd4_cd8_S_INFg_IL2_N_EM_TEMRA.pdf"

# DEPENDENCIES:
require(dplyr)
require(ggplot2)
require(scales)
require(ggbeeswarm)

# <gg4_s_all> all S-reactive IFNg and/or IL2 CD8 phenotypes aggregated
gg8_s_all = ggplot(x_cd8_s, aes(x = factor(vac_num), 
                                y = 100*sum_frequency, 
                                col = COVID_status_per_visit )) + 
  geom_boxplot(width = 1, outlier.size = 0, size = .25)+
  ggbeeswarm::geom_beeswarm(dodge.width=1, size = .1, cex = 1)+
  #geom_point(position = position_dodge(width = .5), size = .25) + 
  scale_y_log10(breaks=c(.01,.1,1,10), labels = c(0.01, 0.1, 1, 10)) + 
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  theme_classic() + 
  annotation_logticks(side = "l", size = .25, 
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) + 
  ylab("")+xlab("") + 
  theme(legend.position = "bottom") + 
  scale_color_manual("",values = c("gray","black")) + 
  #facet_wrap(~factor(vac_num), scale = "free_x") +
  ggpubr::stat_compare_means(label = "p.format",  method = "wilcox.test", size = 2) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  ylab(expression("% CD8 IFN"~ gamma ~"or IL-2")) + 
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 9)) + 
  theme(axis.text.x = element_text()) +
  theme(axis.text.y = element_text( size = 8)) +
  theme(axis.title.y = element_text( size = 8)) +
  coord_cartesian(ylim = c(5E-3, 50)) + 
  theme(legend.position = "none")

# <gg4_s_all> all S-reactive IFNg and/or IL2 CD4 phenotypes aggregated
gg4_s_all = ggplot(x_cd4_s, aes(x = factor(vac_num), 
                                y = 100*sum_frequency, 
                                col = COVID_status_per_visit )) + 
  geom_boxplot(width = 1, outlier.size = 0, size = .25)+
  ggbeeswarm::geom_beeswarm(dodge.width=1, size = .1, cex = 1)+
  #geom_point(position = position_dodge(width = .5), size = .25) + 
  scale_y_log10(breaks=c(.01,.1,1,10), labels = c(0.01, 0.1, 1, 10)) + 
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  theme_classic() + 
  annotation_logticks(side = "l", size = .25, 
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) + 
  ylab("")+xlab("") + 
  theme(legend.position = "bottom") + 
  scale_color_manual("",values = c("gray","black")) + 
  #facet_wrap(~factor(vac_num), scale = "free_x") +
  ggpubr::stat_compare_means(label = "p.format",  method = "wilcox.test", size = 2, label.y= 1) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )+
  ylab(expression("% CD4 IFN"~ gamma ~"or IL-2")) + 
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 9)) + 
  theme(axis.text.x = element_text()) +
  theme(axis.text.y = element_text( size = 8)) +
  theme(axis.title.y = element_text( size = 8)) +
  coord_cartesian(ylim = c(5E-3, 50)) + 
  theme(legend.position = "none")




#### 3-MAJOR SUBSETS, BROKEN OUT
# <gg8_s> all S-reactive IFNg and/or IL2 CD8 phenotypes aggregated by CM, EM, TEMRA subsets
gg8_s = ggplot(x_cd8_s_subset %>% filter(coarse_fct != "Naive"), aes(x = factor(vac_num), 
                                                                     y = 100*sum_frequency, 
                                                                     col = COVID_status_per_visit )) + 
  geom_boxplot(width = 1, outlier.size = 0, size = .25)+
  ggbeeswarm::geom_quasirandom(dodge.width=1, size = .1, cex = .25, width = .2)+
  #geom_point(position = position_dodge(width = .75), size = .25) + 
  scale_y_log10(breaks=c(0.001,.01,.1,1,10), labels = c(0.001,0.01, 0.1, 1, 10)) + 
  theme_classic() + 
  annotation_logticks(side = "l", size = .25, 
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) + 
  ylab("")+xlab("") + 
  theme(legend.position = "bottom") + 
  scale_color_manual("SARS-CoV-2 Prior Infection",values = c("gray","black")) + 
  theme(axis.text = element_text(face = "bold")) + 
  facet_grid(~coarse_fct, scale = "free") +
  ggpubr::stat_compare_means(label = "p.format",  method = "wilcox.test", size = 2.5) + 
  #theme(
  #  strip.background = element_blank(),
  #strip.text = element_blank()
  #)+
  coord_cartesian(ylim = c(1E-4, 50)) + 
  ylab(expression("% CD8 IFN"~ gamma ~"or IL-2")) + 
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 9)) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")

# <gg4_s> all S-reactive IFNg and/or IL2 CD4 phenotypes aggregated by CM, EM, TEMRA subsets
gg4_s = ggplot(x_cd4_s_subset %>% filter(coarse_fct != "Naive"), aes(x = factor(vac_num), 
                                                                     y = 100*sum_frequency, 
                                                                     col = COVID_status_per_visit )) + 
  geom_boxplot(width = 1, outlier.size = 0, size = .25)+
  ggbeeswarm::geom_quasirandom(dodge.width=1, size = .1, cex = .25, width = .2)+
  #geom_point(position = position_dodge(width = .75), size = .25) + 
  scale_y_log10(breaks=c(0.001,.01,.1,1,10), labels = c(0.001,0.01, 0.1, 1, 10)) + 
  theme_classic() + 
  annotation_logticks(side = "l", size = .25, 
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) + 
  ylab("")+xlab("") + 
  theme(legend.position = "bottom") + 
  scale_color_manual("SARS-CoV-2 Prior Infection",values = c("gray","black")) + 
  theme(axis.text = element_text(face = "bold")) + 
  facet_grid(~coarse_fct, scale = "free") +
  ggpubr::stat_compare_means(label = "p.format",  method = "wilcox.test", size = 2.5) + 
  #theme(
  #  strip.background = element_blank(),
  #strip.text = element_blank()
  #)+
  coord_cartesian(ylim = c(1E-4, 50)) + 
  ylab(expression("% CD4 IFN"~ gamma ~"or IL-2")) + 
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 9)) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")

#pdf("/Users/kmayerbl/active/mcelrath_manucript/mansucript_prep/pdf_v4/parts/fig1_cd4_cd8_S_INFg_IL2.pdf", 
pdf(figure_output_1 , 
    width = 2.5,
    height =5)
gridExtra::grid.arrange(gg4_s_all, gg8_s_all, ncol = 1)
dev.off()

#pdf("/Users/kmayerbl/active/mcelrath_manucript/mansucript_prep/pdf_v4/parts/fig1_cd4_cd8_S_INFg_IL2_N_EM_TEMRA.pdf", 
pdf(figure_output_2, 
    width = 6,
    height =5)
gridExtra::grid.arrange(gg4_s, gg8_s, ncol =1)
dev.off()
