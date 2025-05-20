# Purpose: Plot pathway selection intensity for each CCA subtype and panCCA

library(cancereffectsizeR)
library(data.table)
library(ggplot2)

path_effects = fread('output/pathway/panCCA_redefined_pathway_effects.txt')[order(-selection_intensity)]
path_info = fread('output/pathway/redefined_pathway_info.txt')

# Mark cancer gene pathways with stars
path_info[is_cancer_gene_path == T, path_display_name := sub('\\((P\\d+)\\)$', '(\\1\u2605)', path_display_name)]
path_effects[path_info, path_display_name := path_display_name, on = 'path_id']
path_effects[, variant_name := path_display_name]
gg = plot_effects(path_effects, group_by = 'path_display_name', x_title = 'Pathway cancer effect',
             y_title = 'Pathway', legend.position = c(.7, .35), color_by = "#DB382D",
             legend_size_name = 'Nonsilent substitution frequency') +
  theme(axis.text = element_text(family = 'Courier'))

ggsave('figures/pan_path_effects.pdf', gg, width = 6, height = 6, device = cairo_pdf)
