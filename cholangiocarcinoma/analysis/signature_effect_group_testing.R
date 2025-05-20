# File: signature_effect_group_testing.R
# Author: Derek Song
# Date: May 9, 2025
# Purpose: Calculate significant differences in mutational signature effect shares between subtypes.

# Load in signature attributions
effect_list = readRDS(file = 'output/cca_signature_effects.rds')
# ^ from analysis/signature_effect_analysis.R

# get sample-level effect shares for each subtype.
# these are already calculated as proportions of total sample effect!
# i.e. the sum of IHC_sample_effects across all sigs for each sample is 1
IHC_sample_effects = effect_list$iCCA$effect_shares$by_sample
signature_names = names(IHC_sample_effects)[grep("^SBS", names(IHC_sample_effects))]
IHC_sample_effects$cca_type = 'IHC'

# repeat with PHC
PHC_sample_effects = effect_list$pCCA$effect_shares$by_sample
# identical(order(signature_names), order(names(PHC_sample_effects)[-1]))
# # TRUE
PHC_sample_effects$cca_type = 'PHC'

# merge tables
effect_table = rbind(IHC_sample_effects, PHC_sample_effects)

# Significance testing (Mann-Whitney U)
num_tests = length(signature_names)
ihc_phc_results = sapply(1:num_tests, function(x) {
  sig = signature_names[x]
  wilcox.test(
    effect_table[cca_type == 'IHC'][[sig]],
    effect_table[cca_type == 'PHC'][[sig]],
    alternative = 'two.sided',
    exact = FALSE
  )$p.value
})

names(ihc_phc_results) = signature_names
ihc_phc_results = ihc_phc_results[!is.na(ihc_phc_results)] # when both groups have all 0s, you get NA

# 2) IHC vs DCC
# get sample-level effect shares for DCC
DCC_sample_effects = effect_list$dCCA$effect_shares$by_sample
# identical(order(signature_names), order(names(DCC_sample_effects)[-1]))
# # TRUE
DCC_sample_effects$cca_type = 'DCC'

# merge tables
effect_table = rbind(IHC_sample_effects, DCC_sample_effects)

# Significance testing (Mann-Whitney U)
num_tests = length(signature_names)
ihc_dcc_results = sapply(1:num_tests, function(x) {
  sig = signature_names[x]
  wilcox.test(
    effect_table[cca_type == 'IHC'][[sig]],
    effect_table[cca_type == 'DCC'][[sig]],
    alternative = 'two.sided',
    exact = FALSE
  )$p.value
})

names(ihc_dcc_results) = signature_names
ihc_dcc_results = ihc_dcc_results[!is.na(ihc_dcc_results)] # when both groups have all 0s, you get NA

# 2) PHC vs DCC
# merge tables
effect_table = rbind(PHC_sample_effects, DCC_sample_effects)

# Significance testing (Mann-Whitney U)
num_tests = length(signature_names)
phc_dcc_results = sapply(1:num_tests, function(x) {
  sig = signature_names[x]
  wilcox.test(
    effect_table[cca_type == 'PHC'][[sig]],
    effect_table[cca_type == 'DCC'][[sig]],
    alternative = 'two.sided',
    exact = FALSE
  )$p.value
})

names(phc_dcc_results) = signature_names
phc_dcc_results = phc_dcc_results[!is.na(phc_dcc_results)] # when both groups have all 0s, you get NA

# Helper to compute summary stats
summarize_signature_group = function(dt, group_name, sig) {
  x = dt[[sig]]
  list(
    group = group_name,
    signature = sig,
    min = min(x),
    q25 = quantile(x, 0.25),
    median = median(x),
    q75 = quantile(x, 0.75),
    max = max(x),
    mean = mean(x),
    nonzero_mean = mean(x[x > 0])
  )
}

# Function to process each pairwise comparison
generate_comparison_table = function(grp1_dt, grp2_dt, grp1_name, grp2_name, pvals) {
  sigs = names(pvals)
  out = lapply(sigs, function(sig) {
    if (is.na(pvals[sig])) return(NULL)
    stat1 = summarize_signature_group(grp1_dt, grp1_name, sig)
    stat2 = summarize_signature_group(grp2_dt, grp2_name, sig)
    data.table(
      signature = sig,
      comparison = paste(grp1_name, grp2_name, sep = "-"),
      grp1 = grp1_name,
      grp2 = grp2_name,
      pval = pvals[sig],
      grp1.min = stat1$min,
      grp1.25 = stat1$q25,
      grp1.median = stat1$median,
      grp1.75 = stat1$q75,
      grp1.max = stat1$max,
      grp1.mean = stat1$mean,
      grp1.nonzero_mean = stat1$nonzero_mean,
      grp2.min = stat2$min,
      grp2.25 = stat2$q25,
      grp2.median = stat2$median,
      grp2.75 = stat2$q75,
      grp2.max = stat2$max,
      grp2.mean = stat2$mean,
      grp2.nonzero_mean = stat2$nonzero_mean
    )
  })
  rbindlist(out, fill = TRUE)
}

# Run comparisons
sig_cmp_ihc_phc = generate_comparison_table(IHC_sample_effects, PHC_sample_effects, "IHC", "PHC", ihc_phc_results)
sig_cmp_ihc_dcc = generate_comparison_table(IHC_sample_effects, DCC_sample_effects, "IHC", "DCC", ihc_dcc_results)
sig_cmp_phc_dcc = generate_comparison_table(PHC_sample_effects, DCC_sample_effects, "PHC", "DCC", phc_dcc_results)

# Combine
signature_comparisons = rbindlist(list(sig_cmp_ihc_phc, sig_cmp_ihc_dcc, sig_cmp_phc_dcc), use.names = TRUE)

# Add significance label
signature_comparisons[, signif_label := fifelse(pval < 0.001, '***',
                                                fifelse(pval < 0.01, '**',
                                                        fifelse(pval < 0.05, '*', '')))]

fwrite(signature_comparisons, 'output/signature_effect_comparisons.txt', sep = "\t")
