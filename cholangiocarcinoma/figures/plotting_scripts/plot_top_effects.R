library(cancereffectsizeR)
library(data.table)
library(ggplot2)

# Load helper function to annotate noncoding variants.
helpers = new.env()
source('analysis/helpers.R', local = helpers)

cesa = load_cesa('output/cca_cesa.rds')

# Make top gene effects plot
rec_variants = cesa$variants[maf_prevalence > 1, variant_id]

tmp = variant_counts(cesa, rec_variants, by = 'study')[, -"variant_type"][, .SD, 
                                                                          .SDcols = patterns('variant_id|prevalence')]
melted = melt(tmp, id.vars = 'variant_id')[variable != 'total_prevalence']
multi_study_variants = melted[value > 0, uniqueN(variable), by = 'variant_id'][V1 > 1, variant_id]


effects = cesa$selection$all_effects[multi_study_variants, on = 'variant_id']


to_anno = effects[is.na(gene), variant_id]
annotated = helpers$annotate_noncoding(cesa, to_anno, ces.refset.hg38$transcripts)
annotated[gene == 'intergenic', variant_label := variant_name]
annotated[gene != 'intergenic', variant_label := paste0(variant_name, '\n', bare_anno)]


effects[variant_type == 'aac', variant_label := sub('.* ', '', variant_name)]
effects[annotated, let(gene = i.gene, variant_label = i.variant_label), on = 'variant_id']
effects[cesa$variants, essential_splice := essential_splice, on = 'variant_id']

# Because of CES variant prioritization, everything that didn't get annotated already is splice-disrupting.
stopifnot(effects[variant_type == 'snv' & is.na(variant_label), all(essential_splice)])
effects[variant_type == 'snv' & is.na(variant_label), variant_label := paste0(variant_name, '\n(splice)')]


#top_coding = effects[variant_type == 'aac'][order(-selection_intensity), unique(gene)[1:30]]

effects[cesa$variants, conseq := fcase(aa_ref == aa_alt & essential_splice == FALSE, 'Silent',
                                       essential_splice == TRUE, 'Splice-disrupting',
                                       aa_alt == 'STOP', 'Premature stop codon',
                                       aa_ref != aa_alt, 'Missense'), on = 'variant_id']

fp = copy(effects)
fp[is.na(conseq), conseq := 'Other']
fp[, conseq := factor(conseq, levels = c('Premature stop codon', 'Splice-disrupting', 'Missense', 'Other', 'Silent'))]

# Previous color: "#DB382D"
set.seed(9140) # for reproducible ggrepel
gg = plot_effects(fp, group_by = 'gene', topn = 30, 
             label_individual_variants = 'variant_label', color_by = 'conseq', legend_color_name = 'Consequence',
             viridis_option = 'A', x_title = 'Cancer effect') + 
  scale_fill_viridis_d(name = 'Consequence', option = 'B')

ggsave('figures/top_effects.png', gg, width = 9, height = 6.5)

## No longer need to deal with questionable single-study variants.
# to_check = fp[conseq == 'Silent', variant_id]
# select_variants(cesa, variant_ids = to_check, include_subvariants = T)[variant_type == 'snv', unique(variant_id)]

# 1:146994545_C>A, 1:146994545_C>G, 1:146994545_C>T, 7:45084086_C>T, 15:50492857_C>T, 
# 19:21424723_A>G, 19:21972913_G>A, 19:22758865_T>C




# A claim for results section
fgfr2 = fp[gene == 'FGFR2', variant_id]
with_these_fgfr2 = samples_with(cesa, any_of = fgfr2)
cesa$samples[with_these_fgfr2, table(cca_type, fgfr2_fusion, exclude = NULL)]

#         fgfr2_fusion
# cca_type FALSE <NA>
#      IHC    12    3

