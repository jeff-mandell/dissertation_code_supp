# Purpose: Perform group comparisons on signature attributions between CCA subtypes
# Used by: figures/signature_stability_dotplot.R

# Load in required libraries
library(data.table)

# Load in bootstrapped attribution data (see signature_attribution_analysis.R):
# Contains the proportional attribution each signature contributes
# to each sample's total attribution, averaged across boots (370 samples)
signature_attributions = fread('output/final_unblended_signature_weights.txt')

sig_cols = names(signature_attributions)[names(signature_attributions) %like% 'SBS']

# Exclude signatures with no attributions (that is, artifact signatures and non-CCA signatures) from testing
all_zero_sigs = names(which(signature_attributions[, sapply(.SD, function(x) all(x == 0)), 
                                                   .SDcols = sig_cols]))
sig_cols = c(setdiff(sig_cols, all_zero_sigs), 'apobec', 'treatment')
signature_attributions[, apobec := SBS2 + SBS13]
signature_attributions[, treatment := SBS32 + SBS86 + SBS87 + SBS99]
signature_attributions = signature_attributions[, .SD, .SDcols = c(sig_cols, 'Unique_Patient_Identifier', 'cca_type')]

# Function to run Mann-Whitney U on two groups subsetted from signature_attributions
compare_groups = function(grps) {
  stopifnot(is.list(grps), length(grps) == 2, uniqueN(names(grps)) == 2)
  results = sapply(sig_cols, function(x) {
    wilcox.test(
      grps[[1]][[x]],
      grps[[2]][[x]],
      alternative = 'two.sided',
      exact = FALSE
    )$p.value
  })
  
  quantile_and_mean = function(x) c(quantile(x), c(mean = mean(x), nonzero_mean = mean(x[x > 0])))
  grp1_quantiles = as.data.table(t(grps[[1]][, sapply(.SD, quantile_and_mean), .SDcols = sig_cols]), 
                                 keep.rownames = 'signature')
  grp2_quantiles = as.data.table(t(grps[[2]][, sapply(.SD, quantile_and_mean), .SDcols = sig_cols]),
                                 keep.rownames = 'signature')
  quantile_names = c('0%', '25%', '50%', '75%', '100%', 'mean', 'nonzero_mean')
  setnames(grp1_quantiles, quantile_names, paste0('grp1.', c('min', '25', 'median', '75', 'max', 'mean', 'nonzero_mean')))
  setnames(grp2_quantiles, quantile_names, paste0('grp2.', c('min', '25', 'median', '75', 'max', 'mean', 'nonzero_mean')))
  comparison_name = paste(names(grps), collapse = '-')
  all_quantiles = merge.data.table(grp1_quantiles, grp2_quantiles, by = 'signature')
  output = data.table(comparison = comparison_name, signature = sig_cols, 
                      grp1 = names(grps)[1], grp2 = names(grps)[2], pval = results)
  output = merge.data.table(output, all_quantiles, by = 'signature')
  return(output)
}


IHC = signature_attributions[cca_type == 'IHC']
PHC = signature_attributions[cca_type == 'PHC']
DCC = signature_attributions[cca_type == 'DCC']
EHC = signature_attributions[cca_type %in% c('EHC', 'PHC', 'DCC')]

signature_comparisons = rbindlist(list(compare_groups(list(IHC = IHC, EHC = EHC)),
                                       compare_groups(list(IHC = IHC, PHC = PHC)),
                                       compare_groups(list(IHC = IHC, DCC = DCC)),
                                       compare_groups(list(DCC = DCC, PHC = PHC))))

# Test cited in manuscript
signature_comparisons[signature == 'apobec'][comparison %in% c('IHC-PHC', 'DCC-PHC'), .(comparison, pval)]
# comparison        pval
# <char>       <num>
#   1:    IHC-PHC 0.015128034
# 2:    DCC-PHC 0.001738913

# Merge in signature detection frequency from bootstraps
raw_mp_out = fread('output/bootstrapped_mp_out.txt.gz')
raw_mp_out = raw_mp_out[Unique_Patient_Identifier %in% signature_attributions$Unique_Patient_Identifier]
raw_mp_out = raw_mp_out[, .SD, .SDcols = c('Unique_Patient_Identifier', 'boot',
                                           setdiff(sig_cols, c('treatment', 'apobec')))]
raw_mp_out[signature_attributions, cca_type := cca_type, on = 'Unique_Patient_Identifier']

melted = melt(raw_mp_out, id.vars = c('Unique_Patient_Identifier', 'boot', 'cca_type'), 
            variable.factor = FALSE, variable.name = 'signature')
putative = melted[, .(is_detected = mean(value > 0) >= .5), 
                      by = c('Unique_Patient_Identifier', 'signature')][is_detected == T, -"is_detected"]
putative_apobec = unique(putative[signature %in% c('SBS2', 'SBS13'),
                                  .(Unique_Patient_Identifier, signature = 'apobec')])
putative_treatment = unique(putative[signature %in% c('SBS32', 'SBS86', 'SBS87', 'SBS99'),
                                     .(Unique_Patient_Identifier, signature = 'treatment')])
putative = rbind(putative, putative_apobec, putative_treatment)

detection_wp = melted[putative, on = c('Unique_Patient_Identifier', 'signature')][, .(mean_detection_wp = mean(value > 0)), 
                                                                   by = c('Unique_Patient_Identifier', 'signature')]


melted_attr = melt(signature_attributions, id.vars = c('Unique_Patient_Identifier', 'cca_type'), 
                   variable.factor = FALSE, variable.name = 'signature')
attributions_wp = melted_attr[putative, on = c('Unique_Patient_Identifier', 'signature')]


mean_wp = attributions_wp[cca_type != 'EHC', .(mean_wp = mean(value)), by = c('cca_type', 'signature')]
mean_wp_ehc =  attributions_wp[cca_type != 'IHC', .(mean_wp = mean(value)), by = 'signature']
mean_wp_ehc$cca_type = 'EHC'
mean_wp = rbind(mean_wp, mean_wp_ehc)

signature_names = setdiff(unique(signature_comparisons$signature), c('apobec', 'treatment'))
ihc_phc_dcc_detection = raw_mp_out[cca_type %in% c('IHC', 'PHC', 'DCC'), 
                                   lapply(.SD, function(x) mean(x > 0)), 
                                   .SDcols = signature_names, by = 'cca_type']


ehc_detection = raw_mp_out[cca_type %in% c('PHC', 'DCC', 'EHC'),
                           lapply(.SD, function(x) mean(x > 0)),
                           .SDcols = signature_names]
ehc_detection$cca_type = 'EHC'
detection = melt(rbind(ihc_phc_dcc_detection, ehc_detection),
                 id.vars = 'cca_type', value.factor = FALSE, variable.name = 'signature')

detection_wp[signature_attributions, cca_type := cca_type, on = 'Unique_Patient_Identifier']
ihc_phc_dcc_detection_wp = detection_wp[cca_type != 'EHC', .(detection_wp = mean(mean_detection_wp)), 
                                        by = c('signature', 'cca_type')]
ehc_detection_wp = detection_wp[cca_type != 'IHC', .(detection_wp = mean(mean_detection_wp), cca_type = 'EHC'),
                                by = c('signature')]
merged_detection_wp = rbind(ihc_phc_dcc_detection_wp, ehc_detection_wp)

signature_comparisons[detection, grp1.detection := value, on = c(grp1 = 'cca_type', 'signature')]
signature_comparisons[detection, grp2.detection := value, on = c(grp2 = 'cca_type', 'signature')]

signature_comparisons[merged_detection_wp, grp1.detection_wp := detection_wp, on = c(grp1 = 'cca_type', 'signature')]
signature_comparisons[merged_detection_wp, grp2.detection_wp := detection_wp, on = c(grp2 = 'cca_type', 'signature')]


signature_comparisons[mean_wp, grp1.mean_when_present := mean_wp, on = c(grp1 = 'cca_type', 'signature')]
signature_comparisons[mean_wp, grp2.mean_when_present := mean_wp, on = c(grp2 = 'cca_type', 'signature')]

# Significance labels for plots
signature_comparisons[pval < .05, signif_label := '*']
signature_comparisons[pval < .01, signif_label := '**']
signature_comparisons[pval < .001, signif_label := '***']


# Put in percentage of samples with putative detection
putative[signature_attributions, cca_type := cca_type, on = 'Unique_Patient_Identifier']
sample_counts = signature_attributions[, .(total = .N), by = 'cca_type']
ehc_total = sample_counts[cca_type %in% c('EHC', 'PHC', 'DCC'), sum(total)]
sample_counts[cca_type == 'EHC', total := ehc_total]
putative_counts = putative[, .N, by = c('cca_type', 'signature')]
ehc_putative = putative[cca_type %in% c('EHC', 'PHC', 'DCC'), .(.N, cca_type = 'EHC'), by = 'signature']
putative_counts = rbind(putative_counts[cca_type != 'EHC'], ehc_putative)
putative_counts[sample_counts, prop_present := N/total, on = 'cca_type']

signature_comparisons[putative_counts, grp1.prop_present := prop_present, on = c(grp1 = 'cca_type', 'signature')]
signature_comparisons[putative_counts, grp2.prop_present := prop_present, on = c(grp2 = 'cca_type', 'signature')]

total_samples = uniqueN(signature_attributions$Unique_Patient_Identifier)
overall_presence = putative[, .(pan_prop_present = .N/total_samples), by = 'signature']
signature_comparisons[overall_presence, pan_prop_present := pan_prop_present, on = 'signature']


# For manuscript:
overall_presence[order(-pan_prop_present)][1:6]
# signature pan_prop_present
# <char>            <num>
# 1:      SBS1        0.8467742
# 2:    apobec        0.4811828
# 3: treatment        0.4784946
# 4:      SBS2        0.4381720
# 5:     SBS42        0.3602151
# 6:     SBS87        0.3118280

fwrite(signature_comparisons, 'output/signature_comparisons.txt', sep = "\t")
