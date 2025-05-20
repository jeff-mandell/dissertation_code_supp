library(data.table)
library(survival)
library(gtsummary)
library(gt)

stopifnot(packageVersion('gtsummary') >= as.package_version('2.2'))

sample_key = fread('combined_sample_key.txt')
setnames(sample_key, 'Unique_Patient_Identifier', 'patient_id')

pw_info = fread('output/pathway/redefined_pathway_info.txt')

# We'll mark cancer gene pathways with nice stars. (But you may need to remove if they cause rendering problems.)
pw_info[is_cancer_gene_path == T, path_display_name := sub('\\((P\\d+)\\)$', '(\\1\u2605)', path_display_name)]


path_to_gene = fread('output/pathway/redefined_path_to_gene.txt')
needs_ellipses = path_to_gene[, .N, by = 'path_id'][N > 3, path_id]
pw_info[, display_dots := fcase(path_id %in% needs_ellipses, '…', default = '')]
pw_info[, genes := paste0(top_genes, display_dots)]

pw_info[, table_display_name := paste0(stringr::str_pad(paste0('P', pw_rank),
                                                        pad = ' ', width = 8, side = 'right'),
                                       '(', genes, ')')]
pw_info[is_cancer_gene_path == TRUE, table_display_name := sub('(^P\\d+)', '\\1\u2605', table_display_name)]


run_cox = function(features, path_ids, signif_threshold = .1) {
  stopifnot(is.data.table(features),
            all(path_ids %in% names(features)))
  all_samples = features$patient_id
  path_ids = unique(path_ids)
  melted = melt(features[, .SD, .SDcols = c('patient_id', path_ids)], 
                id.vars = 'patient_id', variable.name = 'path_id')
  
  melted = melted[value != ''][, .SD[1], by = c('patient_id', 'path_id')] # not double-counting indel/SBS in same path
  melted[, value := 1]
  melted$patient_id = factor(melted$patient_id, levels = all_samples)
  fm = dcast(melted, patient_id ~ path_id, value.var = 'value', fill = 0)
  missing_samples = data.table(patient_id = setdiff(all_samples, fm$patient_id))
  missing_samples[, (path_ids) := 0]
  fm = rbind(fm, missing_samples)
  fm[, names(.SD) := lapply(.SD, as.numeric), .SDcols = setdiff(names(fm), 'patient_id')]
  fm[sample_key, let(surv_month = surv_month, surv_status = surv_status, pM = pM), on = 'patient_id']
  
  fm = fm[complete.cases(fm)]
  affected_counts = as.list(fm[, colSums(.SD), .SDcols = c(path_ids, 'pM')])
  cox1 = coxph(Surv(surv_month, surv_status) ~ ., data = fm[, -'patient_id'])
  
  stopifnot(cox1$n == fm[, .N])
  
  cox1_coef = as.data.table(summary(cox1)$coefficients, keep.rownames = 'feature')[order(`Pr(>|z|)`)]
  
  not_signif = cox1_coef[`Pr(>|z|)` > signif_threshold, feature]
  cox1_coef = cox1_coef[`Pr(>|z|)` <= signif_threshold]
  pw_disclaimer = paste0(length(not_signif), ' pathways included in the model with p > ', format(signif_threshold), ' not shown.')
  if(length(not_signif) == 0) {
    pw_disclaimer = ''
  }
  ordered_path_ids = setdiff(cox1_coef$feature, 'pM')
  
  table_labels = c(setNames(pw_info$table_display_name, pw_info$path_id),
                   list(pM = 'Presence of\nmetastatic disease'))
  
  # gene_labels = c(setNames(pw_info$genes, pw_info$path_id),
  #                 list(pM = '—'))
  
  curr_table_labels = table_labels[names(table_labels) %in% cox1_coef$feature]
  
  footer = paste0('N = ', cox1$n, ' patients (', cox1$nevent, ' events). ', pw_disclaimer)
  
  # Credit to https://stackoverflow.com/questions/65665465/grouping-rows-in-gtsummary
  # for how to get "Pathway mutation" as a group label
  cox_table = tbl_regression(cox1, label = curr_table_labels, exponentiate = TRUE, 
                              include = cox1_coef$feature, 
                              pvalue_fun = ~style_pvalue(., digits = 2)) |>
    modify_table_body(
      ~ .x |> 
        mutate(#genes = gene_labels[variable],
               num_affected = affected_counts[variable]) |>
        dplyr::relocate(c(num_affected), .before = estimate) |>
        dplyr::bind_rows(
          tibble::tibble(variable="mut_pw", var_label = "mut_pw", row_type="label",
                         label="Pathway mutation in...")) |> 
        dplyr::arrange(factor(variable, levels=c('pM', 'mut_pw', ordered_path_ids)))
    ) |> 
    add_significance_stars(pattern = '{p.value}{stars}', hide_se = T, hide_ci = F, hide_p = F) |>
    modify_column_indent(columns=label, rows = variable %in% path_ids) |> 
    
    modify_header(label = '**Predictor**',
                  #genes = "**Frequently altered genes**",
                  estimate = '**Hazard ratio**',
                  num_affected = "**Number affected**",
                  p.value = "***P*-value**"
                  ) |>
    #modify_column_alignment(columns = 'genes', align = 'left') |> 
    remove_abbreviation() |> 
    modify_source_note(source_note = footer)
  return(cox_table)
}


features = fread('output/landscape/pan_landscape_plotted_features.txt')
accessible_paths = pw_info[weighted_prop_msk >= .25, path_id]
accessible_pw_cox = run_cox(features, accessible_paths)
gtsave(as_gt(accessible_pw_cox), 'output/survival/cancer_gene_pw_cox.pdf')

top_paths = pw_info[pw_rank <= 15, path_id]
nontarget_samples = sample_key[coverage != 'targeted', patient_id]
nontarget_features = features[patient_id %in% nontarget_samples]
top_pw_cox = run_cox(nontarget_features, top_paths)
gtsave(as_gt(top_pw_cox), 'output/survival/top_pw_cox.pdf')

msk_features = fread('output/landscape/msk_landscape_plotted_features.txt')
msk_cox = run_cox(msk_features, accessible_paths)
gtsave(as_gt(msk_cox), 'output/survival/ihc_msk_cox.pdf')


