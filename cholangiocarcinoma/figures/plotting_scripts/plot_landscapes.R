library(data.table)
library(ggplot2)
library(cowplot)

# sample_key = fread('./combined_sample_key_wClusters.txt')
sample_key = fread('combined_sample_key.txt')
setnames(sample_key, 'Unique_Patient_Identifier', 'patient_id')
all_tmb = fread('output/all_sample_in_coverage_TMB.txt')
setnames(all_tmb, 'Unique_Patient_Identifier', 'patient_id')
pw_info = fread('output/pathway/redefined_pathway_info.txt')


# Cap TMB at 10 (will display as "10+").
tmb_min = .1
tmb_max = 10


sample_key[all_tmb, true_tmb := tmb, on = 'patient_id']
sample_key[, tmb := true_tmb]
sample_key[tmb < tmb_min, tmb := tmb_min]
sample_key[tmb > tmb_max, tmb := tmb_max]

# Recode pM as metastatic disease presence
sample_key[, metastatic_disease := fcase(pM == 1, 'Present', pM == 0, 'Absent', default = 'Unknown')]

sample_key[, fgfr2_fusion := fcase(fgfr2_fusion == 1, 'FGFR2 fusion', default = NA)]
sample_key[, other_fusion := fcase(other_fusion == 1, 'Other fusion', default = NA)]

sample_key[, seq_type := fcase(coverage == 'genome', 'WGS',
                       coverage == 'exome', 'WXS',
                       coverage == 'targeted', 'TGS')]

# Set minimum age (will display as "<30").
sample_key[, true_age := age]
sample_key[age < 30, age := 30]
sample_key = sample_key[, .(patient_id, cca_type, age, tmb, true_tmb, metastatic_disease, kit, seq_type, fgfr2_fusion, 
                            other_fusion, surv_month, surv_status)]
sample_key[, cca_subtype := fcase(cca_type == 'IHC', 'Intrahepatic',
                                  cca_type == 'PHC', 'Perihilar',
                                  cca_type == 'DCC', 'Distal',
                                  cca_type == 'EHC', 'Extrahepatic')]

# For each discrete feature, map colors to (all) data values.
color_maps = list(
  metastatic_disease = c("#ef8a62"  = 'Present', "#99d594" = 'Absent', "gray50" = 'Unknown'),
  kit = c("thistle1" = 'MSK341', "skyblue1" = 'MSK410', "olivedrab2" = 'MSK468'),
  cca_subtype = c('dodgerblue4' = 'Intrahepatic',
               'tan' = 'Extrahepatic',
               'darkolivegreen2' = 'Perihilar',
               'gold1' = 'Distal'),
  
  fusions = c("#005824" = 'FGFR2 fusion',  "darkgoldenrod" = 'Other fusion'),
  pathway = c(
    "red4" = 'Copy gain + SBS/Indel',
    "#ca0020" = 'Copy gain',
    "darkblue" = 'Copy decrease + SBS/Indel',
    "#0571b0" = 'Copy decrease',
    "tan4" = 'SBS + Indel',
    "#009E73" = 'Indel',
    "#E69F00" = 'SBS'
  )
)

row_spec = list(seq_type = list(label = 'Sequencing',
                                    legend_name = 'Sequencing type',
                                    color_map =  c(`gold` = 'WXS', `chartreuse4` = 'WGS', `cadetblue` = 'TGS')),
                    cca_subtype = list(label = 'Subtype',
                                    color_map = color_maps$cca_subtype),
                    kit = list(label = 'Panel',
                               legend_name = 'Sequencing panel',
                               color_map = color_maps$kit),
                    age = list(label = 'Age',
                               legend_name = 'Age at diagnosis',
                               viridis_c_opt = list(option = 'A', direction = -1, breaks = c(30, 45, 60, 75), limits = c(25, NA), 
                                                    labels = c('<30', '45', '60', '75'))),
                    tmb = list(label = 'TMB', legend_name = 'TMB',
                               viridis_c_opt = list(option = 'D', trans = 'log10', direction = -1,
                                                    breaks = c(tmb_min, .3, 1, 3, tmb_max), limits = c(tmb_min, tmb_max + 1), 
                                                    labels = c(paste0('<', tmb_min), '0.3', '0.1', '3', paste0(tmb_max, '+')))),
                    metastatic_disease = list(label = 'Metastatic disease',
                                              color_map = color_maps$metastatic_disease),
                    fgfr2_fusion = list(label = 'FGFR2 fusion', legend_name = 'Gene fusions', color_map = color_maps$fusions),
                    other_fusion = list(label = 'Other fusion', color_map = color_maps$fusions, show_guide = FALSE)
)


plot_landscape = function(features, sample_key, path_to_use = pw_info$path_id, 
                          row_spec = row_spec, include_pathways = TRUE) {
  if(rlang::is_scalar_character(features) && file.exists(features)) {
    features = fread(features)
  } else if (is.data.table(features)) {
    features = copy(features)
  } else {
    stop('features expected to be file path or data.table.')
  }
  stopifnot(rlang::is_bool(include_pathways))
  
  dt = copy(features)
  dt = merge.data.table(dt, sample_key, by = "patient_id", all.x = TRUE)
  
  path_to_use = path_to_use[path_to_use %in% names(features)]
  curr_pw_info = pw_info[path_to_use, on = 'path_id']
  pw_names = curr_pw_info$path_display_name
  pathway_descrip = setNames(pw_names, curr_pw_info$path_id)
  
  col_order = path_to_use
  dt[, (col_order) := lapply(.SD, factor, levels = unname(color_maps$pathway)), .SDcols = col_order]
  
  setorderv(dt, c(col_order, 'age'), order = 1, na.last = T)
  tmp = melt(dt[, .SD, .SDcols = c("patient_id", col_order)], id.vars = 'patient_id')
  tmp = tmp[curr_pw_info$path_id, on = 'variable']
  top_path = tmp[! is.na(value), .(top_path = .SD[1, variable]), by = 'patient_id']
  top_path[, cluster_id := match(top_path, col_order)]
  dt[top_path, cluster_id := cluster_id, on = 'patient_id']
  dt[is.na(cluster_id), cluster_id := max(dt$cluster_id, na.rm = T) + 1]
  dt = dt[patient_id %in% features$patient_id]
  dt[, let(x = 1:.N, y = 1)]
  left_label_max_size =  max(nchar(pw_names))
  
  make_row = function(dt, column, label = column, legend_name = label, viridis_c_opt = list(),
                      color_map = NULL, show_guide = TRUE) {
    dt = copy(dt)
    dt[, curr_col := dt[[column]]]
    
    if(! is.list(viridis_c_opt)) {
      stop('viridis_c_opt expected to be type list.')
    }
    if(length(viridis_c_opt) > 0 && ! is.null(color_map)) {
      stop('Use no more than one of arguments viridis_c_opt, color_map.')
    }
    
    # If there's a custom mapping of data values to colors,
    # replace the data values with the colors (and we'll do scale_fill_identity).
    if(! is.null(color_map)) {
      stopifnot(all(na.omit(dt$curr_col) %in% color_map))
      dt[, curr_col := names(color_map)[match(dt$curr_col, color_map)]]
    }
    
    gg = ggplot(dt[, .(x, y, curr_col, cluster_id)], aes(x, y)) + 
      geom_tile(aes(fill = curr_col), show.legend = TRUE) + 
      scale_x_continuous(expand = expansion(0, 0)) + 
      scale_y_continuous(breaks = 1,
                         labels = stringr::str_pad(label, left_label_max_size, 'left'),
                         expand = expansion(0, 0)) +
      facet_grid(cols = vars(cluster_id), scales = 'free_x', space = 'free') +
      theme_minimal() + 
      theme(plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.line = element_blank(), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_text(family = 'Courier'),
            legend.title = element_text(family = 'Courier', size = rel(.8), face = 'bold'),
            panel.grid = element_blank(),
            strip.text = element_blank(), strip.background = element_blank(), panel.spacing = unit(2, 'pt'))
    guide = guide_none()
    if(length(viridis_c_opt) > 0) {
      if(show_guide == T) {
        if(is.null(viridis_c_opt$direction)) {
          to_reverse = F
        } else if(viridis_c_opt$direction == '-1') {
          to_reverse = T
        }
        guide = guide_colorbar(direction = 'horizontal')
      }
      viridis_c_opt = c(viridis_c_opt, list(name = legend_name, guide = guide))
      gg = gg + do.call(scale_fill_viridis_c, viridis_c_opt) + theme(legend.key.width = unit(.75, 'cm'))
      #gg = gg + scale_fill_gradient(name = label, low = colors[1], high = colors[2], trans = 'identity')
    } else if(! is.null(color_map)) {
      if(show_guide == TRUE) {
        guide = guide_legend(direction = 'vertical')
      }
      gg = gg + scale_fill_identity(guide = guide, 
                                    name = legend_name, labels = color_map, limits = names(color_map), 
                                    breaks = names(color_map)) + theme(legend.key.size = unit(8, 'pt'))
    }
    gg = gg + theme(legend.title.position = 'top', legend.justification = c(0, 1),
                    legend.text = element_text(family = 'Courier'))
    if(show_guide) {
      leg = get_plot_component(gg, 'guide-box', return_all = TRUE)[[1]]
    } else {
      # For spacing
      leg = NULL
    }
    
    gg = gg + theme(legend.position = 'none', panel.background = element_rect(fill = 'grey92', color = NA))
    return(list(gg, leg))
  }
  
  subplots = lapply(names(row_spec), function(x) do.call(make_row, c(list(dt = dt, column = x), row_spec[[x]])))
  
  # Add cluster IDs to the cluster color row
  # subplots[[1]][[1]] = subplots[[1]][[1]] + geom_label(data = dt[, .SD[.N/2], by = 'cluster_id'], 
  #                                                      aes(label = cluster_id), 
  #                                                      label.size = .2, label.padding = unit(.2, 'lines'),
  #                                                      label.r = unit(0, 'pt'), size = 3.75, alpha = .85)
  
  if(include_pathways) {
    all_pway_values = unique(na.omit(unlist(dt[, .SD, .SDcols = names(pathway_descrip)])))
    curr_color_map = rev(color_maps$pathway[color_maps$pathway %in% all_pway_values])
    subplots = c(subplots, lapply(names(pathway_descrip), function(x) make_row(dt, x, label = pathway_descrip[x], 
                                                                               color_map = curr_color_map,
                                                                               show_guide = x == names(pathway_descrip)[1],
                                                                               legend_name = 'Pathway mutation')))
  }

  
  row_grid = cowplot::plot_grid(plotlist = lapply(subplots, '[[', 1), ncol = 1)
  all_legends = lapply(subplots, '[[', 2)
  legend_grid = cowplot::plot_grid(plotlist = all_legends[! sapply(all_legends, is.null)], ncol = 1, rel_heights = 1)
  legend_grid = plot_grid(legend_grid, ggplot() + theme_void(), ncol = 1, rel_heights = c(.95, .05))
  
  full_grid = plot_grid(row_grid, ggplot() + theme_void(), legend_grid, nrow = 1, rel_widths = c(.81, .02, .17), align = 'hv', axis = 't')
  final_grid = plot_grid(full_grid, scale = .9)
  final_grid$data = copy(dt)
  return(final_grid)
}

# For MSK, can use the 19 pathways (see prep_landscape.R) that have >25% of present variants covered in MSK.
ihc_effects = fread('output/pathway/subtype_path_effects_with_diff_testing.txt')
ihc_effects = unique(ihc_effects[cca_type == 'IHC'], by = c('path_id'))[, .(path_id, path_display_name, selection_intensity)]
pw_info[ihc_effects, ihc_si := i.selection_intensity, on = 'path_id']
msk_landscape = plot_landscape(features = 'output/landscape/MSK_path_cna_features.txt', sample_key = sample_key, 
                               path_to_use = pw_info[order(-ihc_si), path_id], 
                               row_spec = row_spec[c('kit', 'age', 'tmb', 'metastatic_disease',
                                                     'fgfr2_fusion', 'other_fusion')])
ggsave('figures/landscape/MSK_landscape.pdf', msk_landscape, width = 12, height = 8)
fwrite(msk_landscape$data, 'output/landscape/msk_landscape_plotted_features.txt', sep = "\t")


# ihc_nontarget_landscape = plot_landscape(features = 'output/landscape/nontarget_ihc_features.txt', sample_key = sample_key,
#                                          row_spec = row_spec[c('age', 'tmb', 'metastatic_disease')],
#                                          path_to_use = pw_info[order(-ihc_si)][, .SD[.I %in% 1:25 | weighted_prop_msk > .25]][, path_id])
# ggsave('figures/IHC_nontarget_landscape.png', ihc_nontarget_landscape, width = 12, height = 8)
# fwrite(ihc_nontarget_landscape$data, 'output/landscape/ihc_nontarget_landscape_plotted_features.txt', sep = "\t")

# ihc_landscape = plot_landscape(features = 'output/landscape/IHC_landscape_features.txt',
#                                sample_key = sample_key,
#                                row_spec = row_spec[c('seq_type', 'age', 'tmb', 'metastatic_disease')],
#                                path_to_use = pw_info[order(-ihc_si)][, .SD[.I %in% 1:25 | weighted_prop_msk > .25]][, path_id])
# ggsave('figures/IHC_landscape.png', ihc_landscape, width = 12, height = 8)
# fwrite(ihc_landscape$data, 'output/landscape/ihc_landscape_plotted_features.txt', sep = "\t")

# For panCCA, we'll use top pathways and the 
msk_pw = names(msk_landscape$data[, .SD, .SDcols = patterns('pw')])


pan_pw = pw_info[pw_rank <= 10 | is_cancer_gene_path][order(pw_rank), path_id]

# Put stars in cancer gene path names. Note that you may need a variety of PDF rendering
# pre
pw_info[is_cancer_gene_path == T, path_display_name := sub('\\((P\\d+)\\)$', '(\\1\u2605)', path_display_name)]
pan_landscape = plot_landscape(features = 'output/landscape/pan_landscape_features.txt',
                               sample_key = sample_key,
                               path_to_use = pan_pw,
                               row_spec = row_spec[c('seq_type', 'cca_subtype', 'age', 'tmb', 'metastatic_disease')])
ggsave('figures/landscape/pan_landscape.pdf', pan_landscape, width = 12, height = 8, device = cairo_pdf)
fwrite(pan_landscape$data, 'output/landscape/pan_landscape_plotted_features.txt', sep = "\t")


