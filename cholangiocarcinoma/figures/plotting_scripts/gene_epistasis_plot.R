# Purpose: Plot epistatic interactions between pairs of genes across CCA types.

library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cancereffectsizeR)

msk_exome_output_panCCA = fread(file='output/epistasis_output/gene_epi_MSKexome_panCCA.csv')
exome_output_panCCA = fread(file='output/epistasis_output/gene_epi_exome_panCCA.csv')

# variant A and B are always alphabetical (variant A < B).
# Excluding the exome inferences for genes with exome + MSK inferences.
all_epistasis = rbindlist(list(exome_tgs = msk_exome_output_panCCA,
                               exome_only = exome_output_panCCA[! msk_exome_output_panCCA, 
                                                                on = c('variant_A', 'variant_B')]), 
                          idcol = 'which_inference')
all_epistasis[, fdr := p.adjust(p_epistasis, method = 'fdr')]
fwrite(all_epistasis, 'output/epistasis_output/panCCA_epistatic_effects.txt', sep = "\t")


# Filter to FDR significance
signif_epistasis = all_epistasis[fdr < .05]
fp = signif_epistasis[nAB/n_total > .02 | AB_epistatic_ratio < 1][order(AB_epistatic_ratio)]


fp[, ab_label := format(round(nAB/expected_nAB_null, 2))]

gg = plot_epistasis(fp, pairs_per_row = 5)[[1]] +
  geom_text(aes(x = x, y = 5e4, label = ab_label), nudge_x = .3)
  
ggsave(file='figures/gene_epistasis_panCCA.png', gg,
       width = 7, height = 5)
