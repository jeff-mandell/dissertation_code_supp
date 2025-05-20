library(data.table)
library(cancereffectsizeR)
library(rtracklayer)

cesa = load_cesa('output/cca_cesa.rds')

# For TMB calculation, we'll use just the variants that are strictly within panel target regions.
maf_info = fread('maf_file_summary.txt')
sk = cesa$samples
sk[maf_info, bed := paste0('targeted_regions/', bed), on = 'kit']

stopifnot(sk[is.na(bed), all(coverage == 'genome')],
          all(file.exists(unique(na.omit(sk$bed)))))

bed_groups = split(sk, sk$bed)
tmb_list = list()
for(bed in names(bed_groups)) {
  input_maf = cesa$maf[Unique_Patient_Identifier %in% bed_groups[[bed]]$Unique_Patient_Identifier]
  curr_bed = import.bed(bed)
  bed_size_mb = sum(width(reduce(unstrand(curr_bed))))/1e6
  maf = preload_maf(maf = input_maf, refset = 'ces.refset.hg38', coverage_intervals_to_check = curr_bed)
  maf = maf[dist_to_coverage_intervals == 0]
  tmb = maf[, .(tmb = .N/bed_size_mb), by = 'Unique_Patient_Identifier']
  tmb_list[[bed]] = tmb
}

all_tmb = rbindlist(tmb_list)

# Get hg38 primary assembly genome size
genome_size = sum(seqlengths(cesa$reference_data$genome)[c(1:22, 'X', 'Y')])


## Note: Jusakul WGS calls consist of SNVs and near-coding indels. (Intergenic indels were removed).
# For comparison with the WXS data sets, we'll calculate TMB over default exome regions.
wg_maf_input = cesa$maf[Unique_Patient_Identifier %in% wg_samples]
wg_maf = preload_maf(maf = wg_maf_input, refset = 'ces.refset.hg38', 
                     coverage_intervals_to_check = ces.refset.hg38$default_exome)
wg_maf = wg_maf[dist_to_coverage_intervals == 0]
default_exome_size_mb = sum(width(reduce(ces.refset.hg38$default_exome)))/1e6
wg_tmb = wg_maf[, .(tmb = .N/default_exome_size_mb), by = 'Unique_Patient_Identifier']
all_tmb = rbind(all_tmb, wg_tmb)

# A few samples are left with a TMB of zero due to all variants lying outside of strict coverage.
remaining_samples = setdiff(sk$Unique_Patient_Identifier, all_tmb$Unique_Patient_Identifier)
stopifnot(length(remaining_samples) < 10)

all_tmb = rbind(all_tmb,
                data.table(Unique_Patient_Identifier = remaining_samples, tmb = 0))
stopifnot(uniqueN(all_tmb) == sk[, .N], uniqueN(all_tmb$Unique_Patient_Identifier) == sk[, .N])

# Annotate where burden just reflects nonsynonymous mutations (the TGS data sets, incidentally).
all_tmb[sk, excludes_synonymous := excludes_synonymous, on = 'Unique_Patient_Identifier']

fwrite(all_tmb, 'output/all_sample_in_coverage_TMB.txt', sep = '\t')


