## Read in MAFs from data sources, lift from hg19 to hg38 where needed,
## apply filters (remove possible germline variants and most variants from repetitive regions),
## and exclude cases of possible sample duplication between data sources. Write out finalized MAFs
## for each study.

library(data.table)
library(ces.refset.hg38)
library(cancereffectsizeR)
# setwd(project directory)

chain_file = 'reference/chains/hg19ToHg38.over.chain' # licensed for non-commercial use
maf_info = fread('maf_file_summary.txt')

maf = rbindlist(lapply(1:nrow(maf_info), function(i) {
  genome_build = maf_info$input_build[i]
  stopifnot(genome_build %in% c('hg19', 'hg38'))
  if(genome_build == 'hg19') {
    curr_chain_file = chain_file
  } else {
    curr_chain_file = NULL
  }
  if(is.na(maf_info$bed[i])) {
    curr_bed = NULL
  } else {
    curr_bed = paste0('targeted_regions/', maf_info$bed[i])
    stopifnot(file.exists(curr_bed))
  }
  curr_maf = maf_info$input_maf[i]
  preloaded = preload_maf(maf = curr_maf, refset = 'ces.refset.hg38', chain_file = curr_chain_file, 
                          coverage_intervals_to_check = curr_bed)
  preloaded = preloaded[is.na(problem)]
}), idcol = 'file', fill = TRUE)
maf[, file := basename(maf_info$input_maf[file])] # Replace MAF file's row number in maf_info with file name

# Verify that the targeted regions seem correct, and also determine whether we should pad intervals
# (allow calls somewhat outside the intervals). It is common to call variants within 100bp of
# intervals, but some studies trimmed data to coverage intervals. We'll do padding of 100bp except
# in cases where it seems like the data providers trimmed.
coverage_stats =  maf[, .(covered = mean(dist_to_coverage_intervals == 0), 
                            within_100 = mean(dist_to_coverage_intervals <= 100),
                            within_1000 = mean(dist_to_coverage_intervals <= 1000)), by = 'file']

coverage_stats[, covered_regions_padding := 100]
coverage_stats[covered > .99, covered_regions_padding := 0]
coverage_stats = coverage_stats[! is.na(covered)] # exclude the WGS data set, which has no associated BED file
fwrite(coverage_stats, 'data_prep/maf_coverage_stats.txt', sep = "\t")

# Add some more info and save updated version of file
maf_info$covered_regions_padding = NULL
maf_info[, input_maf_basename := basename(input_maf)]
maf_info[coverage_stats, covered_regions_padding := covered_regions_padding, on = c(input_maf_basename = 'file')]
maf_info[covered_regions_padding == 0, covered_regions_name := kit]
maf_info[covered_regions_padding > 0, covered_regions_name := paste0(kit, '_ip', covered_regions_padding)]

maf[maf_info, covered_regions_padding := covered_regions_padding, on = c(file = 'input_maf_basename')]

# Prepare names for final MAF files
maf_info$input_maf_basename = NULL
maf_info[, final_maf := paste0('final_mafs/', basename(input_maf))]
maf_info[, final_maf := sub('maf.gz', 'final.maf.gz', final_maf)]
setcolorder(maf_info, c('input_maf', 'final_maf'))
fwrite(maf_info, 'maf_file_summary.txt', sep = "\t", na = "NA", quote = F)

# Dropping out-of-coverage records.
maf = maf[is.na(covered_regions_padding) | dist_to_coverage_intervals <= covered_regions_padding]

# The WGS data set has a huge number of intergenic indels called: so many that indels comprise ~75% of the combined
# data set. Since the project's core analysis is on coding regions and SNVs, we'll remove all indels that are not near
# coding transcripts.
maf = maf[, .(file, Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele, variant_id, variant_type)]
maf = preload_maf(maf, refset = 'ces.refset.hg38', keep_extra_columns = TRUE,
                  coverage_intervals_to_check = makeGRangesFromDataFrame(ces.refset.hg38$transcripts[type %in% c('UTR', 'CDS')]))
maf = maf[variant_type == 'snv' | dist_to_coverage_intervals < 200]

snv_stats = maf[variant_type == 'snv', .(prop_repetitive_region = mean(repetitive_region), prop_poss_germline = mean(germline_variant_site), 
                                        prop_either = mean(germline_variant_site | repetitive_region)), by = 'file'][order(-prop_either)]
fwrite(snv_stats, 'data_prep/snv_call_filter_stats.txt', sep = "\t")
indel_stats = maf[variant_type != 'snv', .(prop_repetitive_region = mean(repetitive_region), prop_poss_germline = mean(germline_variant_site), 
                                         prop_either = mean(germline_variant_site | repetitive_region)), by = 'file'][order(-prop_either)]
fwrite(indel_stats, 'data_prep/indel_call_filter_stats.txt', sep = "\t")

# Remove possible germline variants and calls in repetitive regions except those with COSMIC mutation tier annotations.
maf = maf[germline_variant_site == F & (repetitive_region == FALSE | cosmic_site_tier %in% 1:3)]

# Check for possible sample duplicates
poss_dups = check_sample_overlap(maf_list = split(maf, maf$file))

# Examine pairs with more than one shared mutation.
to_examine = poss_dups[variants_shared > 1, .(source_A, source_B, variants_A, variants_B, variants_shared, greater_overlap, sample_A, sample_B)]

# We'll trust MSK curation (these are all cases of two shared variants between pairs of TGS samples).
to_examine = to_examine[! (source_A %like% 'MSK' & source_B %like% 'MSK' & variants_shared == 2)]

# We see 10 (out of 15) Chan-On samples appear to have be resequenced in the larger, more recent Jusakul study.
chan_on_exclusions = to_examine[source_A %like% 'wgs_jusakul' & source_B %like% 'chan-on', sample_B]
stopifnot(length(chan_on_exclusions) == 10)

## Previous note on Ong data (data no longer included in this project):
## 5 out of 8 Ong samples appear resequenced in Jusakul WGS. An additional Ong sample is likely
## resequenced in Jusakul TGS. Since it's hard to verify sample duplication with TGS data and
## there's so much overlap between these studies, we'll go ahead and exclude all Ong samples from
## analysis.

# gao_P7251T and gao_P1026T appear to be duplicates. We'll exclude gao_P1026T (which has fewer variants).
# jiao_CHOL11 and nepal_Patient5 appear to also be sequenced in the TCGA cohort. We'll exclude the non-TCGA samples.
# nepal_Patient13 and nepal_Patient6 appear to also be sequenced in Jiao. We'll keep samples from the larger Jiao study.

# The rest of possible duplicates are 2-3 shared variants and are likely chance overlap (or
# reflective of shared within-study calling error.)
other_exclusions = c("gao_P1026T", "jiao_CHOL11", "nepal_Patient5", "nepal_Patient13", "nepal_Patient6")
to_exclude = c(chan_on_exclusions, other_exclusions)
stopifnot(all(to_exclude %in% c(to_examine$sample_A, to_examine$sample_B)))

# Exclude the identified samples and record decisions
stopifnot(all(to_exclude %in% maf$Unique_Patient_Identifier))
maf = maf[! Unique_Patient_Identifier %in% to_exclude]
to_examine[, decision := 'keep both']
to_examine[sample_A %in% to_exclude, decision := 'exclude A']
stopifnot(! to_examine[decision == 'exclude A', any(sample_B %in% to_exclude)])
to_examine[sample_B %in% to_exclude, decision := 'exclude B']
setcolorder(to_examine, c('sample_A', 'sample_B'))
fwrite(to_examine, 'data_prep/CCA_possible_sample_duplicates.txt', sep = "\t")

maf = maf[, .(file, Chromosome, Start_Position, Reference_Allele, Tumor_Allele, Unique_Patient_Identifier)]
maf_info[, input_basename := basename(input_maf)]
maf[maf_info, final_file := final_maf, on = c(file = 'input_basename')]

final_mafs = split(maf[, -c("file", "final_file")], maf$final_file)
for (i in 1:length(final_mafs)) {
  fwrite(final_mafs[[i]], file = names(final_mafs)[i], sep = "\t")
}
