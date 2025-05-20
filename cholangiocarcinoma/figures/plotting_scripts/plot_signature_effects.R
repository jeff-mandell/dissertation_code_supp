# File: plot_signature_effects.R
# Author: Derek Song
# Date: May 2, 2025
# Purpose: Visualize CCA mutational attribution and effect shares to COSMIC signatures.

library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg38)
source('figures/plotting_scripts/plot_signature_effects_helper.R') # plotting helper function

# Load in signature attributions
effect_list = readRDS(file = 'output/cca_signature_effects.rds') # from analyis/signature_effect_analysis.R

# reproduce cannataro signature groupings and colors, but with some customization
signature_groupings = list(
  "Deamination with age, clock-like" = "SBS1",
  "Unknown, clock-like" = "SBS5",
  "APOBEC" = c("SBS2", "SBS13"),
  "Defective homologous recombination" = "SBS3",
  "Tobacco" = c("SBS4", "SBS29", "SBS92"),
  "Chemotherapeutic agents" = c("SBS32", "SBS86", "SBS87", "SBS99"),
  "Alcohol-associated" = "SBS16",
  "Aristolochic acid exposure" = c("SBS22a", "SBS22b"),
  "Occupational haloalkane exposure" = "SBS42",
  "Aflatoxin exposure" = "SBS24",
  "Mismatch repair defects" = c("SBS15", "SBS20", "SBS21", "SBS44")
)

# Color mapping
cannataro_colors = c(
  "Deamination with age, clock-like" = "gray40",
  "Unknown, clock-like" = "gray60",
  "APOBEC" = "#7570b3",
  "Defective homologous recombination" = "#e7298a",
  "Tobacco" = "#a6761d",
  "Chemotherapeutic agents" = "#1b9e77",
  "Aflatoxin exposure" = "#579c9a",
  "Aristolochic acid exposure" = "#66a61e" ,
  "Occupational haloalkane exposure" = "#e6ab02",
  "Alcohol-associated" = "#d95f02",
  "Non-actionable and unknown signatures" = "black",
  "Mismatch repair defects" = "#8b324d"
)

# Build data.table
dt = rbindlist(lapply(names(signature_groupings), function(desc) {
  sigs = signature_groupings[[desc]]
  data.table(
    name = sigs,
    short_name = as.character(sub("SBS", "", sigs)),
    description = desc,
    prioritize = FALSE,
    color = cannataro_colors[[desc]]
  )
}))

gg = plot_cca_signature_effects(
  mutational_effects = effect_list,
  signature_groupings = dt,
  num_sig_groups = uniqueN(dt$description),
  other_color = "black"
)

ggsave(file = 'figures/weights_by_subtype_Cannataro.png', plot = gg, width = 2700, height = 1300, units = "px")

# gg = plot_signature_effects(
#   mutational_effects = effect_list
# )
# 
# ggsave(file = 'figures/weights_by_subtype_auto.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# # create manual table of sigs, which includes the top 5 signatures by CEW as well as
# # APOBEC signatures 2 and 13
# 
# tbl = cosmic_signature_info()
# tbl = tbl[short_name %in% c('2', '13', '40a', '5', '87', '42', '1')]
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list,
#   signature_groupings = tbl
# )
# 
# ggsave(file = 'figures/weights_by_subtype_subset.png', plot = gg, width = 2700, height = 1300, units = "px")

# # IHC EHC
# effect_list = list(iCCA = IHC_mut_effects,
#                    eCCA = combined_EHC_mut_effects)
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list,
#   signature_groupings = 'cannataro'
# )
# 
# ggsave(file = 'figures/weights_IHC_EHC_Cannataro.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list
# )
# 
# ggsave(file = 'figures/weights_IHC_EHC_auto.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list,
#   signature_groupings = tbl
# )
# 
# ggsave(file = 'figures/weights_IHC_EHC_subset.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# # PHC DCC
# effect_list = list(pCCA = PHC_mut_effects,
#                    dCCA = DCC_mut_effects)
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list,
#   signature_groupings = 'cannataro'
# )
# 
# ggsave(file = 'figures/weights_PHC_DCC_Cannataro.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list
# )
# 
# ggsave(file = 'figures/weights_PHC_DCC_auto.png', plot = gg, width = 2700, height = 1300, units = "px")
# 
# gg = plot_signature_effects(
#   mutational_effects = effect_list,
#   signature_groupings = tbl
# )
# 
# ggsave(file = 'figures/weights_PHC_DCC_subset.png', plot = gg, width = 2700, height = 1300, units = "px")
