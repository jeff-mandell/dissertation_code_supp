# Purpose: Run pairwise differential selection by CCA type on all pathways with effect estimates.

library(cancereffectsizeR)
library(data.table)

cesa = load_cesa("output/cca_cesa.rds")
pw_effects = fread('output/pathway/panCCA_redefined_pathway_effects.txt')
pw_info = fread('output/pathway/redefined_pathway_info.txt')

pathway_defs = readRDS(file='output/pathway/redefined_pathway_defs.rds')
unique_paths = pw_effects$path_id

## Run differential selection analysis
compare_types = function(unique_paths, cca_type_1, cca_type_2, name1 = cca_type_1, name2 = cca_type_2,
                         cores = 1) {
  stopifnot(rlang::is_scalar_character(name1),
            rlang::is_scalar_character(name2))
  all_output = pbapply::pblapply(unique_paths, function(curr_path_id) {
    
    suppressMessages({
      comp = CompoundVariantSet(cesa=cesa, variant_id = pathway_defs[[curr_path_id]])
      prev = pbapply::pboptions(type = 'none')
      cesa = ces_variant(cesa = cesa, variants = comp, samples = cesa$samples[cca_type %in% cca_type_1], run_name = 'run1')
      cesa = ces_variant(cesa = cesa, variants = comp, samples = cesa$samples[cca_type %in% cca_type_2], run_name = 'run2')
      cesa = ces_variant(cesa = cesa, variants = comp, samples = cesa$samples[cca_type %in% c(cca_type_1, cca_type_2)], run_name = 'combined')
      pbapply::pboptions(type = prev$type)
    })
    
    output = merge.data.table(cesa$selection$run1[, .(path_id = curr_path_id, selection_intensity, loglikelihood, 
                                                      included_with_variant, included_total, ci_low_95, ci_high_95)],
                              cesa$selection$run2[, .(path_id = curr_path_id, selection_intensity, loglikelihood, 
                                                      included_with_variant, included_total, ci_low_95, ci_high_95)],
                              suffixes = c('.grp1', '.grp2'), by = 'path_id')
    combined_run_output = cesa$selection$combined[, .(selection_intensity, loglikelihood, ci_low_95, ci_high_95,
                                                      path_id = curr_path_id)]
    
    
    effect_out = rbindlist(setNames(list(cesa$selection$run1, cesa$selection$run2), 
                                    c(name1, name2)), idcol = 'cca_type')
    effect_out$path_id = curr_path_id
    setcolorder(effect_out, 'path_id')
    effect_out[, variant_name := path_id]
    
    testing = merge.data.table(output, combined_run_output, by = 'path_id')
    testing[, two_group_lik := loglikelihood.grp1 + loglikelihood.grp2]
    testing[, chisquared := -2 * (loglikelihood - two_group_lik)]
    testing[, p := pchisq(chisquared, df = 1, lower.tail = F)]
    setnames(testing, names(testing), sub('\\.grp1$', paste0('.', name1), names(testing)))
    setnames(testing, names(testing), sub('\\.grp2', paste0('.', name2), names(testing)))
    testing = testing[, .(path_id, selection_intensity, loglikelihood, two_group_lik, chisquared, p)]
    return(list(effect_out, testing))
  }, cl = cores)
  
  effect_out = rbindlist(lapply(all_output, '[[', 1))
  testing = rbindlist(lapply(all_output, '[[', 2))
  return(list(effects = effect_out,
              testing = testing))
}

ihc_ehc = compare_types(unique_paths = unique_paths, cca_type_1 = 'IHC', 
                        cca_type_2 = c('PHC', 'DCC', 'EHC'), name2 = 'EHC', cores = 2)

phc_dcc = compare_types(unique_paths = unique_paths, cca_type_1 = 'PHC', cca_type_2 = 'DCC', cores = 2)

ihc_phc = compare_types(unique_paths = unique_paths, cca_type_1 = 'IHC', cca_type_2 = 'PHC', cores = 2)
ihc_dcc = compare_types(unique_paths = unique_paths, cca_type_1 = 'IHC', cca_type_2 = 'DCC', cores = 2)


subtype_path_effects = list(ihc_ehc = ihc_ehc, phc_dcc = phc_dcc, ihc_phc = ihc_phc, ihc_dcc = ihc_dcc)
saveRDS(subtype_path_effects, 'output/pathway/subtype_path_effects.rds')


all_effects = rbindlist(lapply(subtype_path_effects, '[[', 1), idcol = 'pair')
all_testing = rbindlist(lapply(subtype_path_effects, '[[', 2), idcol = 'pair')
all_effects[all_testing, let(p_diff_select = p, fdr = fdr), on = c('pair', 'path_id')]
all_effects[pw_info, path_display_name := path_display_name , on = 'path_id']
setcolorder(all_effects, c('pair', 'path_display_name', 'p_diff_select', 'fdr'))
fwrite(all_effects, 'output/pathway/subtype_path_effects_with_diff_testing.txt', sep = "\t")



