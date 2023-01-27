# 2022_12_16_Figure_3
# 2022_12_16_Figure_3_ics_post1_post2.R

# INPUTS:
d_available1 = readr::read_tsv('data/ics_post_1st_post_2nd/available_samples.tsv')
df = readr::read_tsv('data/ics_post_1st_post_2nd/ics_post_1st_post_2nd.tsv')
hla = readr::read_tsv('data/ics_post_1st_post_2nd/hla.tsv')
# OUTPUT:
# Panels used in Figure 3 
figure_output_1 = "figures/2022_12_16_Figure_3_f1_cd8_s1.pdf"
figure_output_2 = "figures/2022_12_16_Figure_3_f1_cd8_s2.pdf"
figure_output_3 = "figures/2022_12_16_Figure_3_f1_cd4_s1.pdf"
figure_output_4 = "figures/2022_12_16_Figure_3_f1_cd4_s2.pdf"

# FUNCTION TO MAKE FIGURE WITH 2 POST-IMMUNIZATION TIMEPOINT 
# The selection criteria is based on subjects with blood draws between 6-30 days
# after the 1st and 2nd mRNA immunization.
figure_1 <- function(tsub = "CD8+", 
                     antigens =  c("COV2 CON S1"), 
                     cyto = "IFNg+", 
                     ylab_text ='% CD8+ T cells expressing IFNg (S1 pool)\n(background-subtracted)'  ){
  # POINTS IN TIMERANGE FOR SELECTION 1
  sel_2 = d_available1 %>% 
    #filter(ICS_run == "ICS Run") %>%
    filter(dbdv2 > 5, dbdv2 <= 30) %>%
    mutate(v2 = dbdv1 - dbdv2) %>%
    select(vax1, vax2,v2, infect, pubid, visitno,  dbdv1, dbdv2, sym_v1, hla_a1, hla_a2) %>% 
    mutate(visit_code= 2)
  # POINTS IN TIMERANGE FOR SELECTION 1
  sel_1 <- d_available1 %>% 
    #filter(ICS_run == "ICS Run") %>%
    filter(dbdv2 <=0, dbdv1 > 5, dbdv1 <= 30) %>%
    mutate(v2 = dbdv1 - dbdv2) %>%
    select(v2, vax1, vax2, infect, pubid, visitno, dbdv1, dbdv2, sym_v1, hla_a1, hla_a2) %>% 
    mutate(visit_code= 1)
  samples = bind_rows(sel_1,sel_2)
  
  s1 = samples %>% left_join(df, c("pubid", "visitno")) %>% 
    filter(tcellsub == tsub) %>% 
    filter(antigen %in%  antigens) %>% #c("COV2 S1-V", "COV2 CON S1")) %>% 
    filter(cytokine == cyto) %>% 
    filter(infect == "PRIOR COVID-19 DIAGNOSIS") %>% 
    group_by(pubid, visitno, visit_code, dbdv1, dbdv2, v2, hla_a1, hla_a2, vax1, vax2) %>% 
    summarise(pctpos_adj = median(pctpos_adj), n = n(), nant = length(unique(antigen)), response.mimosa = max(response.mimosa)) %>% 
    select(pubid,v2, visitno, n, nant, pctpos_adj, visit_code, response.mimosa, hla_a1, hla_a2,vax1, vax2) %>% ungroup() %>%
    left_join(hla, by = "pubid") %>% 
    mutate(hla_a1 = ifelse(is.na(hla_a1),A_1, hla_a1)) %>% 
    mutate(hla_a2 = ifelse(is.na(hla_a2),A_2, hla_a2)) %>% 
    mutate(A03 = stringr::str_detect(hla_a1, pattern = "A\\*03") | stringr::str_detect(hla_a2, pattern = "A\\*03")) %>%
    mutate(hla = paste0(hla_a1,"\n",hla_a2))
  
  a = s1 %>% filter(visit_code == 1)
  b = s1 %>% filter(visit_code == 2)
  ab = a %>% left_join(b, by = 'pubid') %>%
    filter(!is.na(visitno.y )) %>% 
    filter(!is.na(visitno.x)) 
  
  
  rank_sum_p_value = wilcox.test(ab$pctpos_adj.x, ab$pctpos_adj.y, paired = TRUE)$p.value
  signed_rank_p_value = wilcox.test(ab$pctpos_adj.x, ab$pctpos_adj.y, paired = TRUE)$p.value
  print(ab)
  print(wilcox.test(ab$pctpos_adj.x, ab$pctpos_adj.y, paired = TRUE))
  gg_df = dplyr::bind_rows(a,b) %>% 
    filter(pubid %in% (ab$pubid)) %>%
    mutate(pctpos_adj_plot = ifelse(pctpos_adj < 0.001, 0.001, pctpos_adj)) %>% ungroup()
  
  
  responders = gg_df %>% group_by(visit_code) %>% 
    summarise(n=n(), s = sum(response.mimosa))  %>% 
    mutate(rsp= paste0(s,'/',n))
  rsp1 = responders[[1,'rsp']]
  rsp2 = responders[[2,'rsp']]
  
  gg_data = gg_df %>% 
    ggplot(aes(x = factor(visit_code), y = pctpos_adj_plot, group = pubid)) + 
    geom_point(aes(shape = vax1), size = .7) +
    scale_shape_manual("",values=c('Moderna' = 20,  'Pfizer' = 1 ))+ 
    geom_line(alpha = .9, size = .1)+
    theme_classic() + 
    #geom_boxplot(outlier.size = 0, aes(group = factor(visit_code)), width = .1 ,alpha = .3)+
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_discrete(expand = c(.1, .1))+
    annotation_logticks(side = 'l', size=.1) +
    coord_cartesian(ylim = c(0.001, 100)) + 
    theme(axis.line.x=element_line(size=.3)) + 
    theme(axis.line.y=element_line(size=.3))+
    theme(axis.title.x = element_blank()) + 
    theme(axis.text.x = element_blank()) + 
    theme(axis.title.y = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) + 
    ylab(ylab_text) +
    theme(axis.ticks =element_line(size=.1))+
    theme(legend.position = "none")#legend.text=element_text(size=5)) #+ guides(shape=guide_legend(nrow=2))
  
  if (tsub == "CD8+"){
    
    gg_data = gg_data + 
      #geom_text(data = filter(gg_df, visit_code == 2), aes(label = hla), size= .5, col = "blue", nudge_x = .15)+
      #geom_text(data = filter(gg_df, visit_code == 1), aes(label = pubid), size= .5, col = "blue", nudge_x = -.1)+ 
      geom_point(data = filter(gg_df, visit_code == 2, A03 == TRUE), pch = "-", size= 2, col = 'red' , position = position_nudge(x = .05))
    print(filter(gg_df, visit_code == 2, A03 == TRUE))
  }
  
  
  
  labels = tibble::tribble(~cat, ~pos, ~series,
                           "Prior COVID-19",'+','1',
                           "mRNA Dose 1",'+','1',
                           "mRNA Dose 2",'-','1',
                           "Prior COVID-19",'+','2',
                           "mRNA Dose 1",'+','2',
                           "mRNA Dose 2",'+','2')%>% 
    mutate(cat = factor(cat, levels = c( "mRNA Dose 2",  "mRNA Dose 1","Prior COVID-19")))
  
  gg_bottom = ggplot(labels, aes(x = series, y = cat, shape = pos)) + 
    geom_text(aes(label = pos), size  =3) + 
    theme_minimal() + 
    scale_x_discrete(expand =c(0.1,0.1))+
    theme(panel.grid = element_blank()) + ylab("") + xlab("") + 
    #theme(panel.border = element_rect(colour = "black", fill=NA, size=0)) + 
    theme(axis.title.y = element_blank()) + 
    theme(axis.text.y = element_text(size = 6)) +
    theme(axis.text.x = element_blank()) +
    theme(plot.caption = element_text(size = 5))
  
  
  p_value_1_2 = paste0('p = ', round(signed_rank_p_value,4))
  gg_data2 = gg_data + 
    geom_segment(data = data.frame( x = 1.1, y = 50, xend = 1.9, yend = 50), aes(group = NULL, x =x, y = y, xend = xend, yend = yend), size = .2 ) +
    geom_text(data = data.frame(x = 1.5, y = 70, label = p_value_1_2) , aes(group = NULL, x = x, y = y, label = label) , size =3 ) + 
    geom_text(data = data.frame(x = 1.1, y = 30, label = rsp1) , aes(group = NULL, x = x, y = y, label = label) , size =3 )+
    geom_text(data = data.frame(x = 1.9, y = 30, label = rsp2) , aes(group = NULL, x = x, y = y, label = label) , size = 3)
  lm = matrix(c(1,1,1,1,1,2,2))
  final_plot = gridExtra::grid.arrange(gg_data2 + theme(plot.margin = margin(t = 10, r = 10, b = 0, l = 15, unit = "pt")), 
                                       gg_bottom + theme(plot.margin = margin(t = 0, r = 10, b = 10, l = 5, unit = "pt")) + 
                                         labs(caption = "p-value, Wilcox Signed Rank Test,\nPaired, Two-Sided"), 
                                       layout_matrix = lm)
  
  average <- gg_df %>% group_by(pubid) %>% summarise(ave = mean(pctpos_adj_plot)) 
  print(average)
  print(gg_df %>% left_join(average, by = "pubid"))
  gg_time = 
    gg_df %>% left_join(average, by = "pubid") %>%
    ggplot(aes(x = dbdv2, y = pctpos_adj_plot, group = pubid)) + 
    geom_point(size = .5) +
    #geom_point(aes(x = v2, y = ave), col = 'blue', pch = 1)+
    geom_line(alpha = .9)+
    theme_classic() + 
    #geom_boxplot(outlier.size = 0, aes(group = factor(visit_code)), width = .1 ,alpha = .3)+
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(expand = c(.1, .1), breaks = c(-28,-21,-14,-7,0,7,14,21,28))+
    annotation_logticks(side = 'l') +
    coord_cartesian(ylim = c(0.001, 100), xlim = c(-30, 30)) + 
    theme(axis.text.x = element_text(size = 8)) + 
    theme(axis.title.y = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 8)) +
    ylab(ylab_text) 
  
  return(list(final_plot = final_plot, timeplot = gg_time, gg_df = gg_df))
}

# THESE ARE SET ABOVE, HERE FOR REFERENCE
#figure_output_1 = "figures/2022_12_16_Figure_3_f1_cd8_s1.pdf"
#figure_output_2 = "figures/2022_12_16_Figure_3_f1_cd8_s2.pdf"
#figure_output_3 = "figures/2022_12_16_Figure_3_f1_cd4_s1.pdf"
#figure_output_4 = "figures/2022_12_16_Figure_3_f1_cd4_s2.pdf"

pdf(figure_output_1, width = 2, height = 4)
f1_cd8_s1 = figure_1(tsub = "CD8+", 
                     antigens =  c("COV2 CON S1"),
                     cyto = "IFNg_OR_IL2",
                     ylab_text ='% CD8+ T cells expressing IFNg or IL2 (S1 pool)\n(background-subtracted)'  )
f1_cd8_s1$final_plot
dev.off()

pdf(figure_output_2, width = 2, height = 4)
f1_cd8_s2 = figure_1(tsub = "CD8+",
                     antigens =  c( "COV2 CON S2"),
                     cyto = "IFNg_OR_IL2",
                     ylab_text ='% CD8+ T cells expressing IFNg or IL2 (S2 pool)\n(background-subtracted)'  )
f1_cd8_s2$final_plot
dev.off()

pdf(figure_output_3, width = 2, height = 4)
f1_cd4_s1 = figure_1(tsub = "CD4+",
                     antigens =  c( "COV2 CON S1"),
                     cyto = "IFNg_OR_IL2",
                     ylab_text ='% CD4+ T cells expressing IFNg or IL2 (S1 pool)\n(background-subtracted)'  )
f1_cd4_s1$final_plot
dev.off()

pdf(figure_output_4, width = 2, height = 4)
f1_cd4_s2 = figure_1(tsub = "CD4+",
                     antigens =  c("COV2 CON S2"),
                     cyto = "IFNg_OR_IL2",
                     ylab_text ='% CD4+ T cells expressing IFNg or IL2 (S2 pool)\n(background-subtracted)'  )
f1_cd4_s2$final_plot
dev.off()
