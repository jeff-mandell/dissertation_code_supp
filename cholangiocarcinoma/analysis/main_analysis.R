# Purpose: Perform mutation signature extraction and estimation of neutral gene mutation rates.

# Load required libraries
library(data.table)
library(cancereffectsizeR)
library(ces.refset.hg38)

stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('2.10.2'))
stopifnot(packageVersion('ces.refset.hg38') >= as.package_version('1.3.0'))

# Run from project directory.
# MAFs and BED files are in final_mafs and targeted_regions subdirectories.
files = fread('maf_file_summary.txt')
files[! is.na(bed), bed_path := paste0('targeted_regions/', bed)]
files[, maf_name := gsub('\\.final.*', '', basename(final_maf))]

cesa = CESAnalysis('ces.refset.hg38')
for (i in 1:files[,.N]) {
  coverage = files$coverage[i]
  maf = files$final_maf[i]
  if(coverage == 'genome') {
    cesa = load_maf(cesa, maf = maf, coverage = coverage, maf_name = files$maf_name[i])
  } else {
    cesa = load_maf(cesa, maf = maf, coverage = coverage, 
                    covered_regions = files$bed_path[i],
                    covered_regions_name = files$covered_regions_name[i],
                    covered_regions_padding = files$covered_regions_padding[i],
                    maf_name = files$maf_name[i])
  }
}

# Get relative SNV mutation rates across trinucleotide contexts with a mutational signatures method.
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'Biliary-AdenoCA', 
                                                           treatment_naive = F, quiet = T)
set.seed(1010)
n_boot = 1000
boot_output = list()
for(i in 1:n_boot) {
  boot_output[[i]] = trinuc_mutation_rates(
    cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
    signature_exclusions = signature_exclusions, bootstrap_mutations = TRUE,
  )$mutational_signatures$raw_attributions
  if(i %% 10 == 0) {
    message('Finished ', i, ' bootstraps.')
  }
}
combined_mp_out = rbindlist(boot_output, idcol = 'boot')
fwrite(combined_mp_out, 'output/bootstrapped_mp_out.txt.gz', sep = "\t")

# For each sample, calculate the average attribution of every signature across bootstraps
signature_names = names(combined_mp_out)[names(combined_mp_out) %like% 'SBS']
attributions_avg = combined_mp_out[, lapply(.SD, mean), by = "Unique_Patient_Identifier", .SDcols = signature_names]

# Additionally, calculate the proportion that each signature contributes to the total attribution in
# the sample.
attributions_avg[, total := rowSums(.SD), .SDcols = signature_names]
prop_cols = paste0(signature_names, '_prop')
attributions_avg[, (prop_cols) := lapply(.SD, function(x)
  x / total), .SDcols = signature_names]
attributions_avg$total = NULL
fwrite(attributions_avg, 'output/bootstrapped_avg_attributions.txt', sep = "\t")

# Load the proportional weights into the CESAnalysis. First, remove the non-proportional weight columns and
# rename proportional values to match signature names (e.g., SBS1_prop to SBS1).
weights_for_analysis = attributions_avg[, .SD, .SDcols = setdiff(names(attributions_avg), signature_names)]
setnames(weights_for_analysis, names(weights_for_analysis), sub('_prop', '', names(weights_for_analysis)))

# Get group-average weights for samples with at least 50 non-recurrent mutations.
snv_counts_by_trinuc = trinuc_snv_counts(cesa$maf, exclude_recurrent = TRUE, genome = ces.refset.hg38$genome)
snv_counts = colSums(snv_counts_by_trinuc)
well_mutated_samples = intersect(names(snv_counts[snv_counts > 49]), attributions_avg$Unique_Patient_Identifier)
group_avg_weights = weights_for_analysis[Unique_Patient_Identifier %in% well_mutated_samples, 
                                         lapply(.SD, mean), .SDcols = signature_names]
samples_needing_group_average = setdiff(cesa$samples$Unique_Patient_Identifier, attributions_avg$Unique_Patient_Identifier)

# A few samples have all signature weights going to artifact signatures. These samples also need group-average weights.
artifact_signatures = ces.refset.hg38$signatures$COSMIC_v3.4$meta[Likely_Artifact == T, name]
all_artifact_samples = weights_for_analysis[, .(all_artifact = abs(rowSums(.SD) - 1) < 1e-6),.SDcols = artifact_signatures, 
                                            by = 'Unique_Patient_Identifier'][all_artifact == T, Unique_Patient_Identifier]
samples_needing_group_average = union(samples_needing_group_average, all_artifact_samples)
weights_group_averaged_samples = rbindlist(lapply(samples_needing_group_average, 
                 function(x) cbind(data.table(Unique_Patient_Identifier = x), group_avg_weights)))

weights_for_analysis_final = rbind(weights_for_analysis[! samples_needing_group_average, 
                                                        on = "Unique_Patient_Identifier"], 
                                   weights_group_averaged_samples)

cesa = set_signature_weights(cesa, signature_set = 'COSMIC_v3.4', weights = weights_for_analysis_final)

for_wm_anno = cesa$samples[, .(Unique_Patient_Identifier, 
                               sig_analysis_eligible = Unique_Patient_Identifier %in% well_mutated_samples)]
cesa = load_sample_data(cesa, sample_data = for_wm_anno)

# Load sample-specific annotations
sample_key = rbindlist(sapply(list.files(path = 'study_sample_keys/', pattern = 'sample_key.txt', full.names = TRUE),
                              fread, na.strings = c('NA', ''), USE.NAMES = TRUE), fill = T, idcol = 'study')
sample_key[, study := sub('_sample_key\\.txt$', '', basename(study))]
sample_key[cesa$samples, kit := sub('_ip.*', '', covered_regions), on = 'Unique_Patient_Identifier']

stopifnot(sample_key[, .N == uniqueN(Unique_Patient_Identifier)],
          all(cesa$samples$Unique_Patient_Identifier %in% sample_key$Unique_Patient_Identifier))
sample_key = sample_key[, .(Unique_Patient_Identifier, cca_type, age, sex, fluke_status, pM, fgfr2_fusion, other_fusion, 
                            surv_status, surv_month, excludes_synonymous, kit, study)]
sample_key[is.na(excludes_synonymous), excludes_synonymous := FALSE]
cesa = load_sample_data(cesa, sample_key)

# Save a combined sample key
sample_key = cesa$samples[, -'sig_analysis_grp']
fwrite(sample_key, 'combined_sample_key.txt', sep = "\t", na = 'NA')

# For convenience, save a copy of final signature attributions for just the well-mutated samples.
# Group differences will be significance-tested and put in signature stability dotplot
final_weights_wm = cesa$mutational_signatures$biological_weights[well_mutated_samples, on = 'Unique_Patient_Identifier']
final_weights_wm[cesa$samples, cca_type := cca_type, on = "Unique_Patient_Identifier"]
final_weights_wm = final_weights_wm[, .SD, .SDcols = c('Unique_Patient_Identifier', 'cca_type', signature_names)]
fwrite(final_weights_wm, 'output/final_unblended_signature_weights.txt', sep = "\t")


# Estimate neutral gene mutation rates using general-purpose hg38 covariates provided by
# dNdScv developer Inigo Martincorena.
env_for_load = new.env()
load("reference/covariates_hg19_hg38_epigenome_pcawg.rda", envir = env_for_load)
covariates = env_for_load$covs # .rda file provides covariates in an object called covs
rm(env_for_load)

# Run separately for each subtype
cesa =  gene_mutation_rates(cesa, covariates = covariates, samples = cesa$samples[cca_type == 'IHC', Unique_Patient_Identifier])
cesa =  gene_mutation_rates(cesa, covariates = covariates, samples = cesa$samples[cca_type == 'PHC', Unique_Patient_Identifier])
cesa =  gene_mutation_rates(cesa, covariates = covariates, samples = cesa$samples[cca_type == 'DCC', Unique_Patient_Identifier])

# For unspecified EHC, load in average of PHC, DCC
ehc_rates = cesa$gene_rates[, .(rate = (rate_grp_2 + rate_grp_3)/2), by = 'gene']
cesa = set_gene_rates(cesa, rates = ehc_rates, samples = cesa$samples[cca_type == 'EHC', Unique_Patient_Identifier])


# Adjust cores as needed for your system.
cesa = ces_variant(cesa, variants = cesa$variants, run_name = 'all_effects', cores = 4)
save_cesa(cesa, 'output/cca_cesa.rds')

counts_by_type = variant_counts(cesa, by = 'cca_type')
ihc_variants = counts_by_type[IHC_prevalence > 0, variant_id]
dcc_variants = counts_by_type[DCC_prevalence > 0, variant_id]
phc_variants = counts_by_type[PHC_prevalence > 0, variant_id]

all_variants = cesa$variants
setkey(all_variants, 'variant_id')

cesa = ces_variant(cesa = cesa, variants = all_variants[dcc_variants], samples = cesa$samples[cca_type == 'DCC'], run_name = 'DCC', cores = 4)
cesa = ces_variant(cesa = cesa, variants = all_variants[phc_variants], samples = cesa$samples[cca_type == 'PHC'], run_name = 'PHC', cores = 4)
cesa = ces_variant(cesa = cesa, variants = all_variants[ihc_variants], samples = cesa$samples[cca_type == 'IHC'], run_name = 'IHC', cores = 4)

# Save for use in further analyses
save_cesa(cesa, 'output/cca_cesa.rds')




