# Run from dbs subdirectory of dissertation project repo

library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggstatsplot)
library(irr) # for calculation of intraclass correlation

stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('3.0.0.9000'))

# Get TCGA SKCM MAF data. Naming it to make clear that we're including metastatic samples,
# but excluding patients that have multiple sequenced samples.  
mela_maf = 'data/TCGA-SKCM_essential_columns.txt'
if (! file.exists(mela_maf)) {
  tmp_maf = paste0(tempfile(), '.maf')
  get_TCGA_project_MAF(project = 'TCGA-SKCM', filename = tmp_maf, exclude_TCGA_nonprimary = FALSE)
  
  # 2 patients have 2 samples each. For simplicity, we'll remove these patients.
  maf_to_edit = fread(tmp_maf)
  stopifnot(maf_to_edit[multisample_patient == T, 
                        uniqueN(Tumor_Sample_Barcode) == 4 && uniqueN(patient_id) == 2])
  maf_to_edit = maf_to_edit[multisample_patient == FALSE]
  maf_to_edit = maf_to_edit[, .(patient_id, Chromosome, Start_Position, End_Position, Reference_Allele, 
                                Tumor_Seq_Allele2, patient_id, Tumor_Sample_Barcode)]
  fwrite(maf_to_edit, mela_maf, sep = "\t")
  unlink(tmp_maf)
}


# Load SKCM data into CESAnalysis
mela = preload_maf(maf = mela_maf, refset = 'ces.refset.hg38')
cesa = CESAnalysis('ces.refset.hg38')
cesa = load_maf(cesa = cesa, maf = mela)


## Implement dNdScv's indel model, except just for DBS

# Organize covariates data
gr_genes = ces.refset.hg38$gr_genes.dndscv
cv = ces.refset.hg38$covariates$SKCM$rotation
current_cv_genes = rownames(cv)
present_genes = intersect(gr_genes$names, rownames(cv)) # some may be missing, which gets dealt with later
cv = cv[present_genes, ]

# Load "known_cancergenes" from dNdScv
data(list=sprintf("cancergenes_cgc81"), package="dndscv")
cancer_genes = as.data.table(gr_genes)[known_cancergenes, unique(names), on = 'names', nomatch = NULL]

# Annotate DBS with overlapping protein_id
dbs_gr = makeGRangesFromDataFrame(cesa$maf[variant_type == 'dbs', .(seqnames = Chromosome,
                                                                    start = Start_Position, end = Start_Position + 1, variant_id,
                                                                    patient_id)],
                                  keep.extra.columns = TRUE)
overlaps = findOverlaps(dbs_gr, gr_genes, type = 'within')


dbs_by_gene = cbind(as.data.table(dbs_gr[queryHits(overlaps)]), as.data.table(gr_genes[subjectHits(overlaps)]))

# We'll accept one DBS per patient per gene
dbs_by_gene = unique(dbs_by_gene[, .(gene = names, patient_id)])
dbs_by_gene = dbs_by_gene[, .N, by = 'gene']

# Add in rows for genes that don't have any DBS
zero_dbs_genes = data.table(gene = setdiff(gr_genes$names, dbs_by_gene$gene), N = 0)
dbs_by_gene = rbind(dbs_by_gene, zero_dbs_genes)

refcds = ces.refset.hg38$RefCDS.dndscv

# Assemble data for model: Gene name, number of DBS appearing in the CDS (coding sequence; called n_indused to match
# dNdScv's code), length of the CDS, and covariates, cancer gene status.
for_model = dbs_by_gene[, .(gene, n_indused = N, 
                           cds_length = sapply(refcds[gene], '[[', 'CDS_length'))]
for_model[, is_known_cancer := gene %in% known_cancergenes] # known_cancergenes from dNdScv

unif_site_rate = for_model[is_known_cancer == FALSE, sum(n_indused)/sum(cds_length)]
for_model[, exp_unif := cds_length * unif_site_rate]

## Leave out genes with no covariates.
for_model = for_model[gene %in% rownames(cv)]
covs = cv[for_model$gene,]

# Cancer genes are excluded from model fitting. All genes will then get predictions from the fitted model.
# nbrdf is for "negative binomial regression data frame" (as in dNdScv)
nbrdf_all = cbind(for_model[, .(n_indused, exp_unif)], covs)
nbrdf = nbrdf_all[! for_model$is_known_cancer]

nb = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf)

# average mutation rates per DBS site within each pid
for_model[, all_rates := exp(predict(nb,nbrdf_all))]
num_samples = uniqueN(cesa$maf$patient_id)

# 9 DBS are possible per CDS site (each site can be first base of a DBS, and there are 9 DBS for any dinucleotide).
for_model[, final_rate := all_rates / (cds_length * 9) / num_samples]

# Fill out rates for all protein_id by gene.
for_model[, pid := sapply(refcds[gene], '[[', 'protein_id')]

dbs_rates = rbindlist(lapply(ces.refset.hg38$RefCDS, '[', c('protein_id', 'real_gene_name')))
setnames(dbs_rates, 'real_gene_name', 'gene')
dbs_rates[for_model, dbs_rate := final_rate, on = 'gene']

# For now, going to drop genes with no rates. (Later, could use nearby genes.)
dbs_rates = dbs_rates[! is.na(dbs_rate)]


# Get relative trinuc rates across the cohort. Add pseudocounts because a couple contexts don't appear.
present_dbs = cesa$maf[variant_type == 'dbs'][! duplicated(variant_id)]
present_dbs[cesa$annotations$dbs, cosmic_dbs_class := cosmic_dbs_class, on = c(variant_id = 'dbs_id')]

# We are counting each DBS once (the handful of recurrent DBS are enriched for selection)
# cancereffectsizeR supplies us with the full set of DBS categories, cosmic_dbs_classes
cosmic_class_count = table(factor(present_dbs$cosmic_dbs_class, levels = cosmic_dbs_classes))
if(0 %in% cosmic_class_count) {
  cosmic_class_count = cosmic_class_count + 1
}
dbs_prop = cosmic_class_count / sum(cosmic_class_count)

# Likewise, pulling internal tabulation of dinucleotide contexts per gene from refset.
dbs_exposure = cancereffectsizeR:::get_ref_data(cesa, 'cds_dbs_exposure')

# Get context-specific mutation rates for each gene.
gene_site_rates = mapply(
  function(pid, avg_rate) {
    (avg_rate * dbs_prop) / sum(dbs_prop * dbs_exposure[[pid]])
  },
  dbs_rates$protein_id , dbs_rates$dbs_rate, SIMPLIFY = FALSE)
names(gene_site_rates) = dbs_rates$protein_id

# Likelihood function for clonal selection model
dbs_lik = function(rate, num_with, num_without) {
  fn = function(gamma) {
    gamma = unname(gamma) # math faster on unnamed vectors
    sum_log_lik = 0
    if (num_without > 0) {
      sum_log_lik = -1 * gamma * num_without * rate
    }
    if (num_with > 0) {
      sum_log_lik = sum_log_lik + sum(log(1 - exp(-1 * gamma * rep.int(rate, num_with))))
    }
    
    # convert to negative loglikelihood and return
    return(-1 * sum_log_lik)
  }
  # Set default values for gamma (SI), which ces_variant will use to set starting value of optimization
  formals(fn)[["gamma"]] = 1
  bbmle::parnames(fn) = "selection_intensity"
  return(fn)
}

## Collect DBS codon changes and run through selection inference.
dbs_codon_change_ids = cesa$variants[variant_type == 'dbs_aac', variant_id]
easy_codon_change = cesa$annotations$dbs_aac[dbs_aac_id %in% dbs_codon_change_ids][pid %in% names(gene_site_rates)]
  

# Get specific rates for DBS codon changes (summing multiple DBS where necessary)
dbs_key = copy(cesa$annotations$aac_dbs_key)[dbs_aac_id %in% easy_codon_change$dbs_aac_id]
dbs_key[cesa$annotations$dbs, cosmic_dbs_class := cosmic_dbs_class, on = 'dbs_id']
dbs_key[easy_codon_change, pid := pid, on = 'dbs_aac_id']

dbs_key[, rate := mapply(function(x, y) gene_site_rates[[x]][y], dbs_key$pid, dbs_key$cosmic_dbs_class)]
dbs_count = cesa$maf[variant_type == 'dbs', .N, by = 'variant_id'] 
dbs_key[dbs_count, N := N, on = c(dbs_id = 'variant_id')]
dbs_key[is.na(N), N := 0]
info_by_id = dbs_key[, .(N = sum(N), rate = sum(rate)), by = 'dbs_aac_id']

# Run the selection inference. May take a minute or two.
dbs_res = rbindlist(mapply(
  function(curr_dbs, num_with, only_rate) {
    num_without = num_samples - num_with
    fn = dbs_lik(only_rate, num_with, num_without)
    par_init = formals(fn)[[1]]
    names(par_init) = bbmle::parnames(fn)
    # find optimized selection intensities
    # the selection intensity for any stage that has 0 variants will be on the lower boundary; will muffle the associated warning
    withCallingHandlers(
      {
        fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9)
      },
      warning = function(w) {
        if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
          invokeRestart("muffleWarning")
        }
        if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
          # a little dangerous to muffle, but so far these warnings are
          # quite rare and have been harmless
          invokeRestart("muffleWarning") 
        }
      }
    )
    
    selection_intensity = bbmle::coef(fit)
    loglikelihood = as.numeric(bbmle::logLik(fit))
    ci = univariate_si_conf_ints(fit, fn, min_si = 1e-3, max_si = 1e9, conf = .95)
    return(list(variant_id = curr_dbs, selection_intensity = selection_intensity, ll = loglikelihood,
                num_with = num_with, rate = only_rate, ci_low_95 = ci[[1]], ci_high_95 = ci[[2]]))
    
  }, info_by_id$dbs_aac_id, info_by_id$N, info_by_id$rate, SIMPLIFY = FALSE))

# Annotate output.
dbs_res[easy_codon_change, c("gene", "aachange") := .(gene, aachange), on = c(variant_id = 'dbs_aac_id')]
dbs_res[, sbs_aac := paste0(gene, '_', aachange)] # not all are possible
dbs_res$held_out = 0 # currently required column for plot_effects(); it's also the correct value
dbs_res[, included_with_variant := num_with]
dbs_res[, included_total := num_samples]
dbs_res[, variant_type := 'dbs']
dbs_res[, protein_id := sub('.*_', '', variant_id)]
dbs_res[ces.refset.hg38$transcripts, is_mane := is_mane, on = 'protein_id']
dbs_res[, variant_name := fcase(is_mane == T, paste0(gene, ' ', aachange),
                                  default = variant_id)]
dbs_res = dbs_res[num_with > 0]
fwrite(dbs_res, 'dbs_effect_inferences.txt', sep = "\t")

recurrent_results = dbs_res[num_with > 1]
median_inference = dbs_res[, .(variant_id = 'variant', variant_name = '(median DBS)',
                               held_out = 0, variant_type = 'dbs', included_with_variant = 1, 
                               included_total = num_samples, selection_intensity = median(selection_intensity),
                               ci_low_95 = median(ci_low_95), 
                               ci_high_95 = median(ci_high_95))]
for_plot = rbind(recurrent_results, median_inference, fill = TRUE)
dbs_effect_plot = plot_effects(for_plot, color_by = "#DB382D", topn = for_plot[, .N],
                               y_title = 'Cancer effect for DBS-induced amino acid change',
                               legend.position = c(.80, .22))
  #ggtitle('Cancer effects of double-base amino acid changes in TCGA SKCM')
ggsave(filename = 'dbs_aac_effects.png', dbs_effect_plot, width = 7, height = 5.5)



# See what DBS codon change variants also have SBS AACs in the data set.
# 183 of 8,656 DBS results can be checked against SBS.
can_check = intersect(dbs_res$variant_name, cesa$variants[variant_type == 'aac', variant_name])

# Run steps for SBS selection inference. Only running the variants that we need to compare to DBS.
cesa = gene_mutation_rates(cesa, covariates = 'SKCM')
signature_exclusions = suggest_cosmic_signature_exclusions('SKCM', treatment_naive = T, quiet = T)
cesa = trinuc_mutation_rates(cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4, 
                             signature_exclusions = signature_exclusions)
cesa = ces_variant(cesa, variants = cesa$variants[can_check, on = 'variant_name'])


# Merge SBS and DBS inferences and compare
dbs_subset = dbs_res[can_check, on = 'variant_name']
sbs_vs_dbs = merge.data.table(dbs_subset, cesa$selection$selection.1, by = 'variant_name',
                          suffixes = c('.dbs', '.sbs'))

sbs_vs_dbs = sbs_vs_dbs[, .(variant_name, dbs_id = variant_id.dbs, N_dbs = included_with_variant.dbs, 
                            N_snv = included_with_variant.sbs, si_dbs = selection_intensity.dbs, 
                            si_snv = selection_intensity.sbs)]
sbs_vs_dbs[cesa$variants, is_syn := aa_ref == aa_alt, on = 'variant_name']

# We'll label recurrent DBS and anything getting high enough effect to have space for a label.
sbs_vs_dbs[, color := 'gray70']
sbs_vs_dbs[N_dbs > 1, color := 'black']
sbs_vs_dbs[, label := variant_name]
sbs_vs_dbs[, gene := gsub('_.*', '', variant_name)]
sbs_vs_dbs[si_snv < 1e3 & si_dbs < 1e4 & N_dbs == 1, label := '']
sbs_vs_dbs[label != '', label := gsub('_', ' ', label)]

# Verify numbers cited in dissertation.
stopifnot(sbs_vs_dbs[, .N] == 183,
          all(sbs_vs_dbs$is_syn == FALSE))

label_fn = function(x) format(x, big.mark = ",", scientific = F)
sbs_dbs_scatter = ggscatterstats(data = sbs_vs_dbs, type = 'nonparametric',
                    point.args = list(size = 2.5, stroke = 0, color = sbs_vs_dbs$color), 
                    x = si_snv, y = si_dbs,  marginal = F, bf.message = F) + 
  scale_x_log10(labels = label_fn, breaks = 10^(1:6)) +
  scale_y_log10(labels = label_fn, breaks = 10^(2:6), limits = c(99, NA)) +
  geom_text_repel(aes(label = label), size = 3, nudge_x = .5, nudge_y = .04, seed = 900) +
  xlab('SBS-induced amino acid substitution') + ylab('DBS-induced amino acid substitution') + 
  ggtitle('Effects of amino acid changes by nucleotide substitution type') + 
  theme_classic() + coord_fixed() +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 16))

# Picking and choosing which correlation test outputs to include in subtitle
sbs_dbs_scatter$labels$subtitle = as.call(as.list(sbs_dbs_scatter$labels$subtitle)[c(1, 4:6, 3)])
ggsave(plot = sbs_dbs_scatter, filename = 'sbs_dbs_scatter.png', width = 8, height = 8, units = 'in')



## For comparison, let's look at cases where distinct SBS induce the same codon change as each other.
# Find AAC that have multiple distinct SBS in the data
all_aac = cesa$variants[variant_type == 'aac', unique(variant_id)]
multi_snv_aac = cesa$annotations$aac_sbs_key[all_aac, .(sbs_id, multi = length(sbs_id) > 1), 
                                            on = 'aac_id', by = 'aac_id'][multi == TRUE, .(aac_id, sbs_id)]
multi_snv_aac[, in_maf := sbs_id %in% cesa$maf$variant_id]

# There are no AACs that are inflicted by more than two distinct SBS in the data set.
triple_check = multi_snv_aac[, .(triple = sum(in_maf) > 2), by = 'aac_id']
stopifnot(all(triple_check$triple == F))

to_use = multi_snv_aac[, .(to_use = sum(in_maf) == 2), by = 'aac_id'][to_use == TRUE]
snv_aac_to_use = multi_snv_aac[in_maf == TRUE & aac_id %in% to_use$aac_id]

# Verify that we have gotten exactly 2 SBS per AAC
stopifnot(snv_aac_to_use[, .N, by = 'aac_id'][, unique(N)] == 2)

sbs_aac_variants = select_variants(cesa, variant_ids = snv_aac_to_use$sbs_id)
cesa = ces_variant(cesa, variants = sbs_aac_variants, run_name = 'sbs_aac')
sbs_aac_output = copy(cesa$selection$sbs_aac)
sbs_aac_output[snv_aac_to_use, aac_id := aac_id, on = c(variant_id = 'sbs_id')]

# Verify that every codon (that is, pid/aa_pos) is only used in one pair of variants
sbs_aac_output[cesa$variants, aa_pos := aa_pos, on = c(aac_id = 'variant_id')]
sbs_aac_output[, pid := gsub('_.*', '', aac_id)]
stopifnot(uniqueN(sbs_aac_output[, .(pid, aa_pos)]) == sbs_aac_output[, .N/2])


# For intraclass correlation, doesn't matter how variants are divided into grp1/grp2
grp1 = sbs_aac_output[, .SD[1], by = 'aac_id']
grp2 = sbs_aac_output[, .SD[2], by = 'aac_id']
stopifnot(identical(grp1$aac_id, grp2$aac_id))

# Lots of these "codon changes" are actually synonymous.
# The synonymous cases provide a sort of negative control.
cesa$variants[sbs_aac_output$aac_id, on = 'variant_id'][, table(aa_ref == aa_alt)]

syn = cesa$variants[sbs_aac_output$aac_id, on = 'variant_id'][aa_ref == aa_alt, unique(variant_id)]
just_syn_grp1 = grp1[syn, on = 'aac_id']
just_syn_grp2 = grp2[syn, on = 'aac_id']
not_syn_grp1 = grp1[! syn, on = 'aac_id']
not_syn_grp2 = grp2[! syn, on = 'aac_id']


nonsyn_sbs = icc(data.table(not_syn_grp1$selection_intensity, not_syn_grp2$selection_intensity), model = 'oneway')
syn_sbs = icc(data.table(just_syn_grp1$selection_intensity, just_syn_grp2$selection_intensity), model = 'oneway')
dbs_sbs = icc(data.table(sbs_vs_dbs[, .(si_dbs, si_snv)]), model = 'oneway')

icc_table = rbindlist(lapply(list(`DBS/SBS same codon change` = dbs_sbs, 
                                  `SBS/SBS same codon change` = nonsyn_sbs,
                                  `SBS same codon synonymous` = syn_sbs), 
                             function(x) {
                               x[c('subjects', 'value', 'lbound', 'ubound')]
                               }), idcol = 'group')
icc_table = icc_table[, .(group, `number of sites` = subjects, 
                          `intra-class correlation` = round(value, 2),
                          `95% CI` = paste0('(', round(lbound, 2), ', ', round(ubound, 2), ')'))]
fwrite(icc_table, 'dbs_sbs_icc_for_table1.txt', sep = "\t")

save_cesa(cesa, 'cesa_for_dbs.rds')

# Additional info for dissertation
prop_syn = cesa$variants[variant_type == 'aac', mean(aa_ref == aa_alt)]
stopifnot(abs(prop_syn - .33) < .02) # roughly 1 in 3

prop_dbs_syn = cesa$variants[variant_type == 'dbs_aac', mean(aa_ref == aa_alt)]
stopifnot(all.equal(prop_dbs_syn, .034, tolerance = .02))
