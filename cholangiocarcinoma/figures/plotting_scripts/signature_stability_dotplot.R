# Purpose: Plot mean contribution and detection frequency of mutational signatures in CCA

library(data.table)
library(ggplot2)


for_plot_whole = fread('output/signature_comparisons.txt')
for_plot_whole[, short_signature := sub('SBS', '', signature)]
for_plot_whole = rbind(for_plot_whole[, .(comparison, short_signature, cca_type = grp1, contribution = grp1.mean, 
                              contribution_wp = grp1.mean_when_present, detection = grp1.detection, 
                              detection_wp = grp1.detection_wp, signif_label,
                              prop_detected = grp1.prop_present)],
                 for_plot_whole[, .(comparison, short_signature, cca_type = grp2, contribution = grp2.mean, 
                              contribution_wp = grp2.mean_when_present, detection = grp2.detection, 
                              detection_wp = grp2.detection_wp, signif_label,
                              prop_detected = grp2.prop_present)])

## Should be in figure
unique(for_plot_whole, by = c('cca_type', 'short_signature'))[short_signature %in% c('32', '86', '87', '99'),
                                                              .(prop_detected, contribution_wp), by = c('cca_type', 'short_signature')]


major_sigs = for_plot_whole[, .(highest_contribution = max(contribution, na.rm = T)),
                            by = 'short_signature'][order(-highest_contribution), short_signature][1:15]

for_plot = for_plot_whole[short_signature %in% major_sigs]

# Get sample counts
signature_attributions = fread('output/final_unblended_signature_weights.txt')
sample_counts_by_type = signature_attributions[, .N, by = 'cca_type']
ehc_count = sample_counts_by_type[cca_type != 'IHC', sum(N)]
sample_counts_by_type[cca_type == 'EHC', N := ehc_count]
sample_counts_by_type[, count_label := paste0('(n = ', N, ')')]
for_plot[sample_counts_by_type, count_label := count_label, on = 'cca_type' ]

# Build subtype labels
for_plot[, type_label := fcase(cca_type == 'IHC', paste0('Intrahepatic\n', count_label),
                             cca_type == 'PHC', paste0('Perihilar\n', count_label),
                             cca_type == 'DCC', paste0('Distal\n', count_label),
                             cca_type == 'EHC', paste0('Extrahepatic\n', count_label))]

for_type_labeling = setNames(unique(for_plot$type_label), unique(for_plot$cca_type))
for_plot[, cca_type := factor(cca_type, levels = c('IHC', 'EHC', 'PHC', 'DCC'))]

# Will denote significance for two sets of group comparisons
signif_labels = unique(for_plot[comparison %in% c('IHC-EHC', 'DCC-PHC'), .(comparison, short_signature, signif_label)])
signif_labels[, x := fcase(comparison == 'IHC-EHC', 1.5,
                           comparison == 'DCC-PHC', 3.5)]


# Sort by (mostly) numerical signature
for_plot[, sig_num := as.numeric(gsub('[^0-9]', '', short_signature))] # remove all non-numeric
for_plot[, sig_suffix := sub('^\\d+', '', short_signature)]
for_plot = for_plot[order(-sig_num, -sig_suffix)]
for_plot[, c('signif_label', 'sig_num', 'sig_suffix', 'comparison') := NULL]
for_plot[, short_signature := factor(short_signature, levels = unique(for_plot$short_signature))]


dotplot = ggplot(data = for_plot, 
                 aes( y = short_signature, x = cca_type, fill = detection_wp, size = contribution)) +
  geom_point(shape = 21, aes(size = contribution_wp, alpha = 'low')) +
  geom_point(shape = 21, aes(size = contribution, alpha = 'high')) +
  scale_alpha_manual(limits = c('low', 'high'), values = c(.2, 1), 
                     labels = c('Mean contribution when\nputatively detected', 'Mean contribution')) +
  geom_text(data = signif_labels, aes(label = signif_label, x = x, y = short_signature), 
            nudge_y = -0.15, size = 7, na.rm = TRUE, inherit.aes = FALSE) +
  scale_x_discrete(labels = for_type_labeling) + scale_y_discrete() +
  scale_fill_viridis_c(option = 'C', labels = scales::label_percent(), 
                       name = 'Detection frequency\n(when putatively detected)') +
  scale_size(range = c(1, 8), labels = scales::label_percent()) +
  labs(x = "Subtype", y = "SBS signature") +
  guides(fill = guide_colorbar(order = 1),
         size = guide_legend(title = 'Mutational contribution',
                             reverse = TRUE, order = 2),
         alpha = guide_legend(title = 'Contribution type', reverse = TRUE, order = 3,
                              override.aes = list(fill = 'orange', size = 3))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 9), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linewidth = 0.05, color = "grey"),
        panel.border = element_blank(), axis.ticks.y = element_blank())

# add up mean contribution
prop_represented = unique(for_plot[, .(cca_type, contribution, short_signature)])[, .(total = sum(contribution)), by = 'cca_type']
# cca_type        total
# <fctr>     <num>
# 1:      IHC 0.7941940
# 2:      DCC 0.8509981
# 3:      EHC 0.8217937
# 4:      PHC 0.8171288

# A manuscript claim: 
stopifnot(min(prop_represented$total) > .79)


# save plot
ggsave(file = 'figures/signature_dotplots.png', plot = dotplot, height = 5.5, width = 9, units = 'in')

