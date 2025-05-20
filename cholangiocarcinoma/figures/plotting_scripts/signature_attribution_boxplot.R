# File: signature_attribution_boxplot.R
# Author: Derek Song
# Date: March 29, 2025
# Purpose: Visualize group comparisons of APOBEC signature attributions between CCA subtypes

# load in libraries
library(cancereffectsizeR)
library(data.table)
library(ggplot2)

# Load in bootstrapped attribution data (see main_analysis.R)
mp_out = fread('output/final_unblended_signature_weights.txt')

set.seed(8142023) # for reproducible jitter

plot_data = mp_out[cca_type %in% c("IHC", "PHC", "DCC")] # n = 365
plot_data$APOBEC = (plot_data$SBS2 + plot_data$SBS13)

# remap labels for proper CCA nomenclature
label_map = c("IHC" = "iCCA", "PHC" = "pCCA", "DCC" = "dCCA")
plot_data[, cca_display := factor(label_map[cca_type], levels = c("iCCA", "pCCA", "dCCA"))]

# Create the plot
gg = ggplot(plot_data, aes(x = cca_display, y = APOBEC, fill = cca_display)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = c("iCCA" = "dodgerblue4", "pCCA" = "darkolivegreen2", "dCCA" = "gold1")) +
  theme_classic() +
  labs(x = "CCA subtype", y = "Proportion") +
  theme(legend.position = "none")

ggsave(file = 'figures/APOBEC_signatures.png', plot = gg, width = 2700, height = 1300, units = 'px')

# # 2) combined EHC and IHC
# mp_out_temp = mp_out
# # add column that lumps DCC and PHC into EHC
# mp_out_temp$cca_type_lumped = ifelse(mp_out_temp$cca_type %in% c("EHC", "DCC", "PHC"), "EHC", "IHC")
# mp_out_temp$APOBEC = (mp_out_temp$SBS2 + mp_out_temp$SBS13)
# 
# gg = ggplot(mp_out_temp, aes(x = cca_type_lumped, y = APOBEC)) +
#   geom_boxplot(aes(fill = cca_type_lumped), outlier.shape = NA) +
#   geom_jitter(width = .3,
#               size = 0.5,
#               color = 'grey39') +
#   scale_fill_manual(breaks = c('EHC', 'IHC'),
#                     values = c("#E69F00", "#0072B2")) +
#   # scale_x_discrete not needed
#   # scale_y_continuous not needed
#   xlab("CCA subtype") +
#   ylab("Proportion of SNVs attributed to SBS2 + SBS13 (APOBEC)") +
#   theme_classic() + theme(legend.position = "none") +
#   theme(plot.margin = margin(0.2, 0.2, 1, 0.8, "cm"))
# 
# ggsave(file = 'figures/EHC_IHC_APOBEC.png', plot = gg)