# 2022-12-16
# 2022_12_16_Figure_4_expanded_clones.R 
# SCRIPT FOR FIGURE 4 
# CODE REFACTORED FOR PUBLICATIONS (USES PUBIDs)
# Examine Association Between Expanded Clones and ICS Data
  # INPUTS:
  # + 'data/tcrb_matched_ics/tcrb_cohort_sample_labels.tsv'
  # + 'data/tcrb_matched_ics/tcrb_cohort_metadata.tsv'
  # + 'data/tcrb_matched_ics/summary_info_for_2x_fold_change_FDR_gt_0.05_expanded_clones.tsv'
  # + 'data/tcrb_matched_ics/tcrb_matched_intercellular_cytokine_staining_IFN_andor_IL2.csv'
  # (see internal_code folder, for requesting ICS data and assigning publication ID)
  # RETURNS: 
  # + 'figures/Figure4_Panels_ICS_and_Expanded_Clones.pdf'
# NOTES:
  # Pre-processing computes a combined percent positive value 
  # combined CD4 and CD8 pct positive based on 
  # weighted on the %CD4 and %CD8 in each sample 
# DEPENDENCIES
require(ggpubr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)
# INPUTS: 
output_figure_name = 'figures/Figure4_Panels_ICS_and_Expanded_Clones.pdf'
sample_info <- readr::read_tsv('data/tcrb_matched_ics/tcrb_cohort_sample_labels.tsv')
sample_meta_data <- readr::read_tsv('data/tcrb_matched_ics/tcrb_cohort_metadata.tsv')
expanded_clones_filename= 'data/tcrb_matched_ics/summary_info_for_2x_fold_change_FDR_gt_0.05_expanded_clones.tsv'
# OUTPUTS:
output_figure_name ='figures/2022_12_16_Figure4_Panels_ICS_and_Expanded_Clones.pdf'




# 1. PREPROCESSING 
# <df> is intercellular cytokine staining (ICS) data, data-set is matched to visit and PUBID of TCRb-repertoires

df = readr::read_csv('data/tcrb_matched_ics/tcrb_matched_intercellular_cytokine_staining_IFN_andor_IL2.csv')
request_cytokine = "IFNg_OR_IL2"
request_antigen  = "COV2 S1-V"
request_tcellsubs = c("CD4+","CD8+")
# <dfr> refers to requested functional data for S1
dfr = df %>%
  filter(cytokine == request_cytokine ) %>% 
  filter(antigen == request_antigen) %>%
  filter(tcellsub %in% request_tcellsubs)
# <df8> and <df4> refer to separate CD8 and CD4 ICS data into two data frames
df8 = dfr %>% filter(tcellsub == "CD8+")
df4 = dfr %>% filter(tcellsub == "CD4+")
# left join <df8> and <df4>, for S1 response
df48_s1 = df8 %>% left_join(df4, by= c("pubid", "visitno"),suffix = c("_cd8","_cd4")) %>% 
  mutate(percent_cd8 = nsub_cd8/(nsub_cd8+nsub_cd4)) %>% 
  mutate(percent_cd4 = 1-percent_cd8) %>%
  mutate(pp48 = (pctpos_adj_cd8 * percent_cd8) + (pctpos_adj_cd4 * percent_cd4)) 

# repeat request but for S2 response
request_cytokine = "IFNg_OR_IL2"
request_antigen  = "COV2 S2-V"
request_tcellsubs = c("CD4+","CD8+")
# <dfr> refers to requested functional data for S2
dfr = df %>%
  filter(cytokine == request_cytokine ) %>% 
  filter(antigen == request_antigen) %>%
  filter(tcellsub %in% request_tcellsubs)
# left join <df8> and <df4>, for S2 response
df8 = dfr %>% filter(tcellsub == "CD8+")
df4 = dfr %>% filter(tcellsub == "CD4+")
df48_s2 = df8 %>% left_join(df4, by= c("pubid", "visitno"),suffix = c("_cd8","_cd4")) %>% 
  mutate(percent_cd8 = nsub_cd8/(nsub_cd8+nsub_cd4)) %>% 
  mutate(percent_cd4 = 1-percent_cd8) %>%
  mutate(pp48 = (pctpos_adj_cd8 * percent_cd8) + (pctpos_adj_cd4 * percent_cd4)) 

# Request variant pools
request_cytokine = "IFNg_OR_IL2"
request_antigen  = "COV2 Var S"
request_tcellsubs = c("CD4+","CD8+")
dfr = df %>%
  filter(cytokine == request_cytokine ) %>% 
  filter(antigen == request_antigen) %>%
  filter(tcellsub %in% request_tcellsubs)
df8 = dfr %>% filter(tcellsub == "CD8+")
df4 = dfr %>% filter(tcellsub == "CD4+")
df48_var = df8 %>% left_join(df4, by= c("pubid", "visitno"),suffix = c("_cd8","_cd4")) %>% 
  mutate(percent_cd8 = nsub_cd8/(nsub_cd8+nsub_cd4)) %>% 
  mutate(percent_cd4 = 1-percent_cd8) %>%
  mutate(pp48 = (pctpos_adj_cd8 * percent_cd8) + (pctpos_adj_cd4 * percent_cd4)) 


# <exp_df>
exp_df = readr::read_tsv(expanded_clones_filename) 
exp_df = readr::read_tsv(expanded_clones_filename) %>% 
  dplyr::filter(type %in% c('post2_pre','post1_pre')) %>% 
  select(post, type, sum_pf_expanded_clones, breadth_expanded_clones, simpson_expanded_clones, shannon_expanded_clones) %>% 
  mutate(breadth_expanded_clones = ifelse(breadth_expanded_clones == 0, 1, breadth_expanded_clones ))
# Left join s1,s2, variant pools by left joining and then combining 
# percent positive CD4,CD8 s1, s2, variants (i.e., pp48_s1_s2_var)
# also divide breadth of exapanded clones by number of productive clones in that sample
df48_s1_s2 = df48_s1 %>% 
  left_join(df48_s2, by = c("pubid","visitno"),  suffix = c("_s1","_s2")) %>% 
  left_join(df48_var, by = c("pubid","visitno"),  suffix = c("","_var"))%>%
  mutate(post = paste0(pubid, "_v",visitno )) %>% 
  mutate(pp48_s1= ifelse(is.na(pp48_s1), 0, pp48_s1 )) %>% 
  mutate(pp48_s2=ifelse(is.na(pp48_s2), 0, pp48_s2 )) %>%    
  mutate(pp48_var=ifelse(is.na(pp48), 0, pp48 )) %>%  
  mutate(pp48_s1_s2 = pp48_s1 + pp48_s2) %>%
  mutate(pp48_s1_s2_var = pp48_s1 + pp48_s2+ pp48_var)%>%
  left_join(exp_df, by = "post") %>% 
  left_join(sample_info, by = "post") %>% 
  mutate(breadth_measure = breadth_expanded_clones / number_of_productive_clones)

df48_s1_s2$number_of_productive_clones

df48_s1_s2$type_num = c("post1_pre" = '1', "post2_pre" = '2')[df48_s1_s2$type]
# Save the checkpoint file
df48_s1_s2 %>% readr::write_tsv('data/tcrb_matched_ics/tcrb_matched_ICS_IFN_andor_IL2_merged_with_expanded_clones.tsv')

# 2. VISUALIZATIONS

# Reload <df48_s1_s2> from checkpoint file 
df48_s1_s2 <- readr::read_tsv('data/tcrb_matched_ics/tcrb_matched_ICS_IFN_andor_IL2_merged_with_expanded_clones.tsv')

# Specify some visual parameters.
colores = c("black","black")
line_col = "black"
line_alpha = .5
line_size = .25
pt_size = .6
add_ptid_labels = T

# Define a theme
theme_comparison = 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = "none") + 
  theme(axis.text = element_text(face = "bold", size = 8)) + 
  theme(axis.title = element_text(face = "bold", size = 8))  

# <gg_func> ggplot graphic showing response to both S1 and S2 
# (combined CD4 and CD8)
gg_func = df48_s1_s2 %>% 
  mutate(sum_pf_expanded_clones = ifelse(sum_pf_expanded_clones< 0.0005, 
                                         0.0005, 
                                         sum_pf_expanded_clones)) %>% 
  ggplot(aes(x =type_num, y = pp48_s1_s2_var)) + 
  geom_line(aes(group = pubid), 
            alpha = line_alpha ,
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type, shape = type, fill = type), size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.1,0.5))+
  xlab("")+ 
  ylab("% IFNg+ or IL2 CD4 and CD8 T cells\nS1 and S2 peptide pools") +
  annotation_logticks(side = 'l', size = .2) + 
  theme_comparison


# <gg_func_S1> ggplot graphic showing functional response to both S1 
# (combined CD4 and CD8)
gg_func_S1 = df48_s1_s2 %>% mutate(sum_pf_expanded_clones = ifelse(sum_pf_expanded_clones< 0.0005, 
                                                                   0.0005, 
                                                                   sum_pf_expanded_clones)) %>% 
  ggplot(aes(x =type_num, y = pp48_s1)) + 
  geom_line(aes(group = pubid),
            alpha = line_alpha, 
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type,
                 shape = type, 
                 fill = type), 
             size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.1,0.5))+
  xlab("")+ 
  ylab("% IFNg+ or IL2 CD4 and CD8 T cells\nS1 peptide pool") +
  annotation_logticks(side = 'l', size = .2) + 
  theme_comparison


# <gg_func_S2> ggplot graphic showing functional response to both S2 
# (combined CD4 and CD8)
gg_func_S2 = df48_s1_s2 %>% 
  mutate(sum_pf_expanded_clones = ifelse(sum_pf_expanded_clones< 0.0005, 
                                         0.0005, 
                                         sum_pf_expanded_clones)) %>% 
  filter(pubid != 'CO001-535869') %>% # CD4 was not measured because of sample limitation
  ggplot(aes(x =type_num, y = pp48_s2)) + 
  geom_line(aes(group = pubid),
            alpha = line_alpha,
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type,
                 shape = type,
                 fill = type),
             size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.1,0.5))+
  xlab("")+ 
  ylab("% IFNg+ or IL2 CD4 and CD8 T cells\nS2 peptide pool") +
  annotation_logticks(side = 'l', size = .2) + 
  theme_comparison

# <gg_freq> ggplot graphic showing frequency of expanded clones
gg_freq = df48_s1_s2 %>% 
  mutate(sum_pf_expanded_clones = ifelse(sum_pf_expanded_clones< 0.0005, 
                                         0.0005, 
                                         sum_pf_expanded_clones)) %>% 
  ggplot(aes(x =type_num, y = sum_pf_expanded_clones*100)) + 
  geom_line(aes(group = pubid),
            alpha = line_alpha,
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type,
                 shape = type,
                 fill = type),
             size = pt_size)+ 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() +
  scale_x_discrete(expand = c(0.1,0.5))+
  xlab("")+ 
  ylab("% TCRb which expanded vs. Pre-Vaccine") +
  annotation_logticks(side = 'l', size = .2) + 
  theme_comparison

# <gg_breadth > ggplot graphic showing breadth of expanded clones
gg_breadth = df48_s1_s2 %>% 
  ggplot(aes(x =type_num, y = breadth_measure)) + 
  geom_line(aes(group = pubid), 
            alpha = line_alpha, 
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type, shape = type, fill = type), size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_x_discrete(expand = c(0.1,0.5))+
  scale_y_log10() + 
  xlab("")+ 
  ylab("Breadth of TCRb Expanded vs. Pre-Vaccine\n(FDR < 0.05)") +
  annotation_logticks(side = 'l', size = .2) + 
  theme_comparison


# <ggS1S2> cross-plot of breadth (x) vs. ICS (y)
ggS1S2_breadth = df48_s1_s2 %>% 
  ggplot(aes(x =100*breadth_measure, 
             y = pp48_s1_s2_var)) + 
  geom_line(aes(group = pubid),
            alpha = line_alpha,
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type,
                 shape = type,
                 fill = type),
             size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() + 
  scale_x_log10() + 
  theme_bw() + 
  annotation_logticks(side = 'lb', size = .2) +  
  #coord_cartesian(xlim = c(0.0001, .35), ylim = c(0.001, 1.5)) + 
  theme(legend.position = "bottom") + 
  xlab("Breadth of TCRb clones expanded vs. Pre-Vaccine") + 
  ylab("% IFNg+ or IL2+ CD4 and CD8 T cells\nS1 and S2 peptide pools") + 
  ggpubr::stat_cor( method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, label.x = -3, label.y =1 , size = 3) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(size = .25)) + 
  theme(axis.text = element_text(face = "bold", size = 8)) + 
  theme(axis.title = element_text(face = "bold", size = 8)) + 
  theme(legend.position = "none")+
  ggrepel::geom_text_repel(data = df48_s1_s2 %>% filter(type == "post2_pre"), aes(label = pubid),
                           force_pull   = 0, # do not pull toward data points
                           nudge_x      = .1,
                           direction    = "y",
                           size = 1.5,
                           angle        = 0,
                           hjust        = 0,
                           segment.size = 0.2,
                           col = "steelblue",
                           max.iter = 1e4, 
                           max.time = 1)

# <ggS1S2_freq> cross-plot of sum of productive frequency (x) vs. ICS (y)
ggS1S2_freq = df48_s1_s2 %>% 
  mutate(sum_pf_expanded_clones = ifelse(sum_pf_expanded_clones< 0.0005,
                                         0.0005, 
                                         sum_pf_expanded_clones)) %>%
  ggplot(aes(x =sum_pf_expanded_clones*100, 
             y = pp48_s1_s2_var)) +
  geom_abline(col = "gray", linetype = "dotted")+
  geom_line(aes(group = pubid),
            alpha = line_alpha,
            size = line_size,
            col = line_col)+
  geom_point(aes(col = type,
                 shape = type,
                 fill = type),
             size = pt_size) + 
  scale_color_manual("",values=colores)+
  scale_fill_manual("",values=colores)+
  scale_shape_manual("",values=c(1,21)) +
  scale_y_log10() + 
  scale_x_log10() + 
  theme_bw() + 
  annotation_logticks(side = 'lb', size = .2) + 
  theme(legend.position = "bottom") + 
  xlab("% TCRb expanded vs. Pre-Vaccine") + 
  ylab("% IFNg+ or IL2+ CD4 and CD8 T cells\nS1 and S2 peptide pools") + 
  ggpubr::stat_cor( method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, label.x = -2, label.y =1 , size = 3) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(size = .25)) + 
  theme(axis.text = element_text(face = "bold", size = 8))+ 
  theme(axis.title = element_text(face = "bold", size = 8)) + 
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(data = df48_s1_s2 %>% filter(type == "post2_pre"), aes(label = pubid),
                           force_pull   = 0, # do not pull toward data points
                           nudge_x      = .1,
                           direction    = "y",
                           size = 1.5,
                           angle        = 0,
                           hjust        = 0,
                           segment.size = 0.2,
                           col = "steelblue",
                           max.iter = 1e4, 
                           max.time = 1)
# PLOT
pdf(output_figure_name,  width = 8, height= 8)
ggarrange(
  ggarrange(gg_func + coord_cartesian(ylim = c(0.01,12)),
            gg_freq +  coord_cartesian(ylim = c(0.01,12)), 
            gg_breadth, 
            ncol = 3, labels = c("A","B","C","D","E")), 
  ggarrange(ggS1S2_freq+ coord_cartesian(ylim = c(0.01,12), xlim =  c(0.01,12)),
            ggS1S2_breadth+ coord_cartesian(ylim = c(0.01,12)),
            ncol = 2, labels = c("F","G")), 
  nrow = 2
)
dev.off()




