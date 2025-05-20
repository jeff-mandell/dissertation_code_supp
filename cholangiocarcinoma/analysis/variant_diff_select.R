library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg38)

cesa = load_cesa("./output/cca_cesa.rds")

# Only need to test variants with prevalence > 2.
variants_to_use = cesa$variants[maf_prevalence > 2, variant_id]


# Test differential selection. Excluding variants with low coverage (<50 in either group) or low prevalence (< 3 combined)
prev_threshold = 3
type_cov_threshold = 20
compare_types = function(variants, cca_type_1, cca_type_2, name1 = cca_type_1, name2 = cca_type_2) {
  
  stopifnot(rlang::is_scalar_character(name1),
            rlang::is_scalar_character(name2))
  
  vc = variant_counts(cesa, variant_ids = variants, by = 'cca_type')
  melted = melt(vc[, -"variant_type"], id.vars = 'variant_id')
  melted[, cca_type := sub('_.*', '', variable)]
  melted = melted[cca_type %in% c(cca_type_1, cca_type_2)]
  prevalences = melted[variable %like% '_prevalence']
  good_prev = prevalences[, sum(value), by = 'variant_id'][V1 >= prev_threshold, variant_id]
  covs = melted[variable %like% '_covering']
  good_cov_1 = covs[cca_type %in% cca_type_1, sum(value), by = 'variant_id'][V1 > type_cov_threshold, variant_id]
  good_cov_2 = covs[cca_type %in% cca_type_2, sum(value), by = 'variant_id'][V1 > type_cov_threshold, variant_id]
  good_cov = intersect(good_cov_1, good_cov_2)
  variants = cesa$variants[variant_id %in% intersect(good_cov, good_prev)]
  stopifnot(variants[, .N] > 0)
  
  cesa = ces_variant(cesa = cesa, variants = variants, samples = cesa$samples[cca_type %in% cca_type_1], run_name = 'run1')
  cesa = ces_variant(cesa = cesa, variants = variants, samples = cesa$samples[cca_type %in% cca_type_2], run_name = 'run2')
  cesa = ces_variant(cesa = cesa, variants = variants, samples = cesa$samples[cca_type %in% c(cca_type_1, cca_type_2)], run_name = 'combined')
  output = merge.data.table(cesa$selection$run1[, .(variant_id, selection_intensity, loglikelihood, 
                                                    included_with_variant, included_total, ci_low_95, ci_high_95)],
                            cesa$selection$run2[, .(variant_id, selection_intensity, loglikelihood, included_with_variant, 
                                                    included_total, ci_low_95, ci_high_95)],
                            suffixes = c('.grp1', '.grp2'), by = 'variant_id')
  combined_run_output = cesa$selection$combined[, .(variant_id, selection_intensity, loglikelihood, ci_low_95, ci_high_95)]
  all_output = merge.data.table(output, combined_run_output, by = 'variant_id')
  all_output[, two_group_lik := loglikelihood.grp1 + loglikelihood.grp2]
  all_output[, chisquared := -2 * (loglikelihood - two_group_lik)]
  all_output[, p := pchisq(chisquared, df = 1, lower.tail = F)]
  
  setnames(all_output, names(all_output), sub('\\.grp1$', paste0('.', name1), names(all_output)))
  setnames(all_output, names(all_output), sub('\\.grp2', paste0('.', name2), names(all_output)))
  effect_out = rbindlist(setNames(list(cesa$selection$run1, cesa$selection$run2), c(name1, name2)), idcol = 'cca_type')
  return(list(effects = effect_out,
              testing = all_output[, .(variant_id, selection_intensity, loglikelihood, two_group_lik, chisquared, p)]))
}

ihc_ehc = compare_types(variants = variants_to_use, cca_type_1 = 'IHC', cca_type_2 = c('PHC', 'DCC', 'EHC'), name2 = 'EHC')
phc_dcc = compare_types(variants_to_use, cca_type_1 = 'PHC', cca_type_2 = 'DCC')
ihc_phc = compare_types(variants_to_use, cca_type_1 = 'IHC', cca_type_2 = 'PHC')
ihc_dcc = compare_types(variants_to_use, cca_type_1 = 'IHC', cca_type_2 = 'DCC')

subtype_variant_effects = list(ihc_ehc = ihc_ehc, phc_dcc = phc_dcc, ihc_phc = ihc_phc, ihc_dcc = ihc_dcc)
saveRDS(subtype_variant_effects, 'output/subtype_variant_effects.rds')


all_effects = rbindlist(lapply(subtype_variant_effects, '[[', 1), idcol = 'pair')
all_testing = rbindlist(lapply(subtype_variant_effects, '[[', 2), idcol = 'pair')
all_effects[all_testing, let(p_diff_select = p, fdr = fdr), on = c('pair', 'variant_id')]
setcolorder(all_effects, c('pair', 'variant_name', 'p_diff_select', 'fdr'))
fwrite(all_effects, 'output/subtype_variant_effects_with_diff_testing.txt', sep = "\t")

