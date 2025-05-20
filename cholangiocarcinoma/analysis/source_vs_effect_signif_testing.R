library(data.table)

mut_effects = readRDS('output/cca_signature_effects.rds')
attributions = fread('output/final_unblended_signature_weights.txt')

all_signatures = setdiff(names(mut_effects$iCCA$effect_shares$by_sample), 'Unique_Patient_Identifier')

source_and_effect = melt(attributions, id.vars = c('Unique_Patient_Identifier', 'cca_type'), 
                         measure.vars = all_signatures, variable.name = 'signature', 
                         variable.factor = FALSE, value.name = 'source_share')

sample_effect_shares = rbindlist(lapply(mut_effects, function(x) x$effect_shares$by_sample), idcol = 'cca_type')
for_merge = melt(sample_effect_shares,  id.vars = c('Unique_Patient_Identifier', 'cca_type'), 
                 measure.vars = all_signatures, variable.name = 'signature', 
                 variable.factor = FALSE, value.name = 'effect_share')
for_merge[, cca_type := fcase(cca_type == 'iCCA', 'IHC',
                              cca_type == 'pCCA', 'PHC',
                              cca_type == 'dCCA', 'DCC', default = NA)]
stopifnot(! anyNA(for_merge$cca_type))

# Note that unspecified EHC samples will have NA effect_share. We're not using these saamples since
# the effect shares are for iCCA, pCCA, dCCA.
source_and_effect[for_merge, effect_share := effect_share, on = c('Unique_Patient_Identifier', 'cca_type')]

# Here's one test
treatment_signatures = c('SBS32', 'SBS86', 'SBS87', 'SBS99')
for_trt = source_and_effect[signature %in% treatment_signatures, 
                                 .(cca_type, 
                                   source_share = sum(source_share), 
                                   effect_share = sum(effect_share),
                                   signature = 'treatment'), by = 'Unique_Patient_Identifier']
apobec_signatures = c('SBS2', 'SBS13')
for_apobec = source_and_effect[signature %in% apobec_signatures, 
                               .(cca_type, 
                                 source_share = sum(source_share), 
                                 effect_share = sum(effect_share),
                                 signature = 'apobec'), by = 'Unique_Patient_Identifier']

for_testing = rbindlist(list(source_and_effect, for_trt, for_apobec), use.names = TRUE)

subtype_testing = for_testing[cca_type != 'EHC', .(pval = wilcox.test(source_share, effect_share, paired = T, exact = F)$p.value),
                             by = c('cca_type', 'signature')]

pan_testing = for_testing[, .(cca_type = 'panCCA', pval = wilcox.test(source_share, effect_share, paired = T, exact = F)$p.value),
                              by = 'signature']
signif_testing = rbind(subtype_testing, pan_testing)
sorted_sbs = c('apobec', 'treatment', all_signatures)
signif_testing = signif_testing[sorted_sbs, on = 'signature']
signif_testing[, fdr := p.adjust(pval, 'fdr')]


# For manuscript
signif_testing[signature %in% c('treatment', 'apobec')]
# cca_type signature         pval          fdr
# <char>    <char>        <num>        <num>
# 1:      IHC    apobec 3.515272e-05 7.030545e-05
# 2:      PHC    apobec 2.287730e-04 4.102136e-04
# 3:      DCC    apobec 1.997860e-10 6.631194e-10
# 4:   panCCA    apobec 2.822131e-14 1.294860e-13
# 5:      IHC treatment 1.395485e-31 2.176957e-30
# 6:      PHC treatment 3.271545e-07 9.279291e-07
# 7:      DCC treatment 1.881012e-13 8.151051e-13
# 8:   panCCA treatment 8.493256e-49 6.624739e-47

fwrite(signif_testing, 'output/source-vs-effect_signif.txt')


