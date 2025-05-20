# File: signature_effect_analysis.R
# Author: Derek Song
# Date: May 2, 2025
# Purpose: Calculate CCA mutational attribution and effect shares to COSMIC signatures.

# load in libraries
library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg38)

# load in CESA with cancer effects calculated (produced from analysis/main_analysis.R)
cesa = load_cesa('output/cca_cesa.rds')

# N = 372 signature analysis eligible samples
eligible_samples = cesa$samples[sig_analysis_eligible == TRUE, Unique_Patient_Identifier]

# calculate source shares and effect shares of mutational signatures.

# iCCA
iCCA_mut_effects = mutational_signature_effects(cesa, effects = cesa$selection$IHC, 
                                                samples = intersect(eligible_samples, cesa$samples[cca_type == 'IHC', Unique_Patient_Identifier]))

# dCCA
dCCA_mut_effects = mutational_signature_effects(cesa, effects = cesa$selection$DCC, 
                                                samples = intersect(eligible_samples, cesa$samples[cca_type == 'DCC', Unique_Patient_Identifier]))

# pCCA
pCCA_mut_effects = mutational_signature_effects(cesa, effects = cesa$selection$PHC, 
                                                samples = intersect(eligible_samples, cesa$samples[cca_type == 'PHC', Unique_Patient_Identifier]))

# # eCCA combined
# combined_EHC_mut_effects = mutational_signature_effects(
#   cesa = cesa,
#   effects = cesa$selection$combined_EHC,
#   samples = mp_out[cca_type %in% c('EHC', 'DCC', 'PHC')]$Unique_Patient_Identifier # 121 samples
# )
# 
# # dCCA + pCCA
# combined_dCCA_pCCA_mut_effects = mutational_signature_effects(
#   cesa = cesa,
#   effects = cesa$selection$combined_DCC_PHC,
#   samples = mp_out[cca_type %in% c('DCC', 'PHC')]$Unique_Patient_Identifier # 111 samples
# )

# Gather all of the mutational_signature_effects() outputs to graph
effect_list = list(iCCA = iCCA_mut_effects,
                   pCCA = pCCA_mut_effects,
                   dCCA = dCCA_mut_effects)

saveRDS(effect_list, file = 'output/cca_signature_effects.rds')

