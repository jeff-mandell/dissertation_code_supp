# Purpose: Run CES on pathways defined by the variants within pathway and consensus exome-intersected genomic ranges.

library(cancereffectsizeR)
library(ces.refset.hg38)
library(data.table)
library(pbapply)

# Adjust as necessary/desired for parallel operations.
cores = 2

stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('2.10.0'))
stopifnot(packageVersion('ces.refset.hg38') >= as.package_version('1.3.0'))

cesa = load_cesa("output/cca_cesa.rds")

# Intersect consensus exome-covered regions with pathway definitions to get
# the regions of each pathway that should be covered in all exome/genome-sequenced samples.
cons_exome = Reduce(intersect, cesa$coverage_ranges$exome) # 23Mb
path_gene_ranges = as.data.table(get_PathScore_coding_regions(genome="hg38"))


cons_exome = Reduce(intersect, cesa$coverage_ranges$exome)
cons_msk = Reduce(intersect, cesa$coverage_ranges$targeted[! names(cesa$coverage_ranges$targeted) == 'Jusakul_TGS'])
msk_exome_inter = intersect(cons_msk, cons_exome) # 0.7 Mb

# For each pathway, merge across genes to get pathway-wide ranges
path_gene_link = fread(file='reference/pathway_info/pathway_gene_link.txt') # path_id, entrez_id
path_ranges = unique(merge.data.table(path_gene_ranges, path_gene_link[, .(path_id, entrez_id)], by = 'entrez_id', all = F, 
                                         allow.cartesian = TRUE)[, .(path_id, seqnames, start, end)])
path_grl = reduce(makeGRangesListFromDataFrame(path_ranges, split.field = 'path_id'))

# Intersect the pathway ranges with consensus exome definitions to get within-coverage pathway ranges.
# May take a couple of minutes; adjust cores as necessary for your machine.
path_grl_exome = pblapply(path_grl, intersect, cons_exome, ignore.strand = T, cl = cores)
names(path_grl_exome) = paste0('pw.', names(path_grl_exome))


path_grl_mskex = pblapply(path_grl, intersect, msk_exome_inter, ignore.strand = T, cl = cores)
names(path_grl_mskex) = paste0('pw.', names(path_grl_mskex))

# Access to elements is faster when not compressed
path_grl_exome = GRangesList(path_grl_exome, compress = F)
path_grl_mskex = GRangesList(path_grl_mskex, compress = F)

# Save a copy (Fun fact: xz is good for GRanges.)
saveRDS(path_grl_exome, 'analysis/path_analysis/pathway_consensus_exome_granges.rds', 
        compress = 'xz')
saveRDS(path_grl_mskex, 'analysis/path_analysis/pathway_consensus_msk-exome_granges.rds', 
        compress = 'xz')

# Get covered sizes of pathways in bases
path_sizes_exome = sum(width(path_grl_exome))
path_sizes_mskex = sum(width(path_grl_mskex))

# Some pathways have no exome coverage in the exome ranges; remove these.
path_grl_exome = path_grl_exome[path_sizes_exome > 0]
path_grl_mskex = path_grl_mskex[path_sizes_exome > 0]
prop_msk = path_sizes_mskex[names(path_grl_mskex)]/path_sizes_exome[names(path_grl_exome)]


# Effects will be estimated across nonsilent variants within consensus-coverage pathway region.
# Subset MAF to nonsilent variants.
maf = cesa$maf

# A little ugly due to an issue where essential_splice could be FALSE or NA (which are equivalent here).
maf[cesa$variants, is_silent := aa_ref == aa_alt & (is.na(essential_splice) | essential_splice == FALSE), 
    on = c(top_consequence = 'variant_name')]

# Only need one record per variant_id
nonsilent_snv_maf = maf[is_silent == FALSE | is.na(is_silent)][variant_type == 'snv'][! duplicated(variant_id)]

snv_gr = makeGRangesFromDataFrame(nonsilent_snv_maf[, .(chr = Chromosome, start = Start_Position, 
                                                        end = Start_Position, variant_id)],
                                  keep.extra.columns = TRUE)

present_variants_by_path = pblapply(path_grl_exome, function(x) snv_gr$variant_id[overlapsAny(snv_gr, x)],
                                    cl = cores)
present_variants_mskex = pblapply(path_grl_mskex, function(x) snv_gr$variant_id[overlapsAny(snv_gr, x)],
                                  cl = cores)


saveRDS(present_variants_by_path, 'output/pathway/initial_exome_pw_defs.rds')
saveRDS(present_variants_mskex, 'output/pathway/initial_tgs_pw_dfs.rds')

# We can use just the MSK definitions for the (few) pathways with prop_msk > .5.
# Then, we can use the regular definitons for the remaining pathways, and we don't lose any variants at all.
which_high = which(prop_msk > .5)
stopifnot(length(which_high) == 237) # 237 out of 3,379
t1 = unname(unlist(present_variants_mskex[which_high]))
not_high_msk = setdiff(1:length(present_variants_by_path), which_high)
t2 = unname(unlist(present_variants_by_path[not_high_msk]))
stopifnot(uniqueN(c(t1, t2)) == uniqueN(unlist(present_variants_by_path)))

# How many genes do we get with our MSK-intersect pathways?
to_check = unique(unname(unlist(present_variants_mskex[which_high])))
gene_check = select_variants(cesa, variant_ids = to_check)
stopifnot(uniqueN(gene_check$gene) == 218)

msk_pw = present_variants_mskex[names(which_high)]
rest_of_pw = present_variants_by_path[! names(present_variants_by_path) %in% names(which_high)]
names(msk_pw) = paste0(names(msk_pw), '.tgs')

initial_pathway_defs = c(msk_pw, rest_of_pw)




## Algorithm to get independent pathway effects.
# Tentative pathways = (List of pathway definitions)
# Final pathways = (Empty list)
# while TRUE:
#   Filter tentative pathways to those mutated in at least 2% of samples.
# Define N = number of tentative pathways.
# Calculate panCCA effects and rank by effect.
# Add pathway 1 to the list of final pathways and remove it from tentative pathways.
# Redefine tentative pathways 2-N with all variants in pathway 1 excluded.

calc_path_effects = function(current_pathway_defs, samples = character(), cores = 1) {
  output = rbindlist(pblapply(names(current_pathway_defs),
                    function(pathway) {
                      prev_opt = pboptions(type = 'none') # locally disable progress bar
                      comp = suppressMessages(CompoundVariantSet(cesa=cesa, variant_id = current_pathway_defs[[pathway]]))
                      pboptions(prev_opt)
                      names(comp) = pathway
                      return(ces_variant(cesa, variants = comp, samples = samples, 
                                           run_name = pathway)$selection[[pathway]])
                    }, cl = cores))
  return(output)
}

freq_threshold = .05
redefine_pathways = function(pathway_defs, samples_to_use, freq_threshold, cores = 1, stop_at = Inf) {
  final_pathway_defs = list()
  current_pathway_defs = copy(pathway_defs)
  final_effects = data.table()
  cached_results = data.table()
  round = 1
  total_pathways = length(current_pathway_defs)
  needing_mut_freq_check = current_pathway_defs
  num_samples = length(samples_to_use)
  
  while(TRUE) {
    prev_opt = pboptions(type = 'none')
    bad_paths = names(needing_mut_freq_check[sapply(needing_mut_freq_check, length) == 0])
    needing_mut_freq_check = needing_mut_freq_check[! names(needing_mut_freq_check) %in% bad_paths]
    
    path_mut_freq = pbsapply(needing_mut_freq_check, 
                             function(x) length(intersect(samples_with(cesa, x), samples_to_use))/num_samples, cl = cores)
    pboptions(prev_opt)
    bad_paths = c(bad_paths, names(path_mut_freq[path_mut_freq < freq_threshold]))
    num_bad = length(bad_paths)
    if(num_bad > 0) {
      message("Removed ", num_bad, " pathways due to low frequency of mutation.")
      total_pathways = total_pathways - num_bad
    }
    current_pathway_defs = current_pathway_defs[! names(current_pathway_defs) %in% bad_paths]
    if(length(current_pathway_defs) == 0) break
    to_run = setdiff(names(current_pathway_defs), cached_results$variant_name)
    if(length(to_run) == 0) {
      curr_effects = cached_results
    } else {
      num_final = final_effects[, .N]
      message("Round ", round, ": ", num_final, " of ", total_pathways, " pathways done. Calculating new effects for ", length(to_run), ' pathways...')
      curr_effects = suppressMessages(calc_path_effects(current_pathway_defs = current_pathway_defs[to_run], samples = samples_to_use, cores = cores))
      curr_effects = rbind(curr_effects, cached_results)
      round = round + 1
    }
    top_path = curr_effects[which.max(selection_intensity)]
    final_effects = rbind(final_effects, top_path)
    top_path_id = top_path$variant_name
    variants_in_top = current_pathway_defs[[top_path_id]]
    final_pathway_defs[[top_path_id]] = variants_in_top
    if(length(final_pathway_defs) >= stop_at) break
    current_pathway_defs[[top_path_id]] = NULL
    previous_sizes = sapply(current_pathway_defs, length)
    current_pathway_defs = lapply(current_pathway_defs, function(x) setdiff(x, variants_in_top))
    new_sizes = sapply(current_pathway_defs, length)
    unchanged = names(current_pathway_defs[previous_sizes == new_sizes])
    needing_mut_freq_check = current_pathway_defs[! names(current_pathway_defs) %in% unchanged]
    
    # Any pathways not changing do not need to be recalculated on the next run
    cached_results = curr_effects[variant_name %in% unchanged]
  }
  return(list(defs = final_pathway_defs, effects = final_effects))
}

tgs_pw_samples = cesa$samples[covered_regions != 'Jusakul_TGS', Unique_Patient_Identifier]
tgs_pw_redefined = redefine_pathways(pathway_defs = msk_pw, samples_to_use = tgs_pw_samples, freq_threshold = .05, cores = 2)

already_used = unique(unname(unlist(tgs_pw_redefined$defs)))
exome_pathway_defs = lapply(rest_of_pw, function(x) setdiff(x, already_used))
nontarget_samples = cesa$samples[coverage != 'targeted', Unique_Patient_Identifier]
exome_pw_redefined = redefine_pathways(pathway_defs = exome_pathway_defs, samples_to_use = nontarget_samples,
                                       freq_threshold = .05, cores = 2, stop_at = 100)

final_pathway_defs = c(tgs_pw_redefined$defs, exome_pw_redefined$defs)
final_effects = rbind(tgs_pw_redefined$effects, exome_pw_redefined$effects)
setnames(final_effects, 'variant_name', 'path_id')

# Save now
saveRDS(list(defs = final_pathway_defs, effects = final_effects), file = 'output/pathway/redefined_pathways.rds')

# Get annotations for the pathway variants
final_path_to_variants = rbindlist(lapply(final_pathway_defs, as.data.table), idcol = 'path_id')
setnames(final_path_to_variants, 'V1', 'variant_id')
stopifnot(uniqueN(final_path_to_variants$variant_id) == final_path_to_variants[, .N])
anno = select_variants(cesa, variant_ids = final_path_to_variants$variant_id, 
                       include_subvariants = T)[variant_type == 'snv']
final_path_to_variants[anno, let(maf_prevalence = maf_prevalence, gene = gene),
                       on = 'variant_id']

variants_in_msk_inter = select_variants(cesa, gr = msk_exome_inter, 
                                        include_subvariants = TRUE)[variant_type == 'snv', variant_id]
final_path_to_variants[, in_msk := variant_id %in% variants_in_msk_inter]
weighted_prop_msk = final_path_to_variants[, .(prop = .SD[in_msk == TRUE, sum(maf_prevalence)]/sum(maf_prevalence)), by = 'path_id']


redefined_path_to_gene = final_path_to_variants[, .(count = sum(maf_prevalence)), by = c('gene', 'path_id')]
redefined_path_to_gene = redefined_path_to_gene[order(-count)]

custom_pw_info = redefined_path_to_gene[, .(gene_list = .(na.omit(gene[1:3]))), by = 'path_id']
custom_pw_info[, top_genes := sapply(gene_list, paste, collapse = ', ')]
custom_pw_info[, gene_list := NULL]

custom_pw_info[final_effects, si := selection_intensity, on = 'path_id']
custom_pw_info[weighted_prop_msk, weighted_prop_msk := prop, on = 'path_id']


pathway_info = fread('reference/pathway_info/pathways_info.txt')
pathway_info[, path_id_for_merge := paste0('pw.', PATH_ID)]
custom_pw_info[, path_id_for_merge := sub('\\.tgs', '', path_id)]
custom_pw_info[pathway_info, let(pathway_name = pathway_name, name_systematic = name_systematic, 
                                 descrip = description_brief), on = 'path_id_for_merge']
custom_pw_info$path_id_for_merge = NULL
custom_pw_info = custom_pw_info[order(-si)]
custom_pw_info[, pw_rank := 1:.N]
custom_pw_info[, path_display_name := paste0(gsub(', ', '|', top_genes), ' (P', 1:.N, ')')]
custom_pw_info[, is_cancer_gene_path := weighted_prop_msk >= .25]

fwrite(custom_pw_info, 'output/pathway/redefined_pathway_info.txt', sep = "\t")
fwrite(redefined_path_to_gene[, .(path_id, gene)], 'output/pathway/redefined_path_to_gene.txt', sep = "\t")
fwrite(final_effects[order(-selection_intensity)], 'output/pathway/panCCA_redefined_pathway_effects.txt', sep = "\t")
saveRDS(final_pathway_defs, 'output/pathway/redefined_pathway_defs.rds')

