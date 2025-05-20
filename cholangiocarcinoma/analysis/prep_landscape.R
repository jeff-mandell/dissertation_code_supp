library(cancereffectsizeR)
library(data.table)

cesa = load_cesa("output/cca_cesa.rds")
present_variants_by_path = readRDS("output/pathway/redefined_pathway_defs.rds")

# Get tables that associate pathway IDs to gene IDs
pathway_info = fread('output/pathway/redefined_pathway_info.txt')
pway_by_gene = fread('output/pathway/redefined_path_to_gene.txt')


prep_landscape_data = function(included_samples, pw_to_use, gene_copy_calls = NULL) {
  maf = cesa$maf[Unique_Patient_Identifier %in% included_samples]
  
  ## Filter to nonsilent variants.
  maf[cesa$variants, is_silent := aa_ref == aa_alt & (is.na(essential_splice) | essential_splice == FALSE), 
      on = c(top_consequence = 'variant_name')]
  
  # Non-SNVs have is_silent = NA.
  nonsilent_maf = maf[is_silent == FALSE | is.na(is_silent)]
  nonsilent_snv_maf = nonsilent_maf[variant_type == 'snv']
  nonsilent_indel_maf = nonsilent_maf[variant_type != 'snv']
  
  suppressMessages({
    gene_mut_by_sample = make_PathScore_input(nonsilent_maf)
    gene_snv_mut_by_sample = make_PathScore_input(nonsilent_snv_maf)
    gene_indel_mut_by_sample = make_PathScore_input(nonsilent_indel_maf)
  })

  gene_mut_by_sample[gene_snv_mut_by_sample, has_snv := TRUE, on = c('patient_id', 'hugo_symbol')]
  gene_mut_by_sample[gene_indel_mut_by_sample, has_indel := TRUE, on = c('patient_id', 'hugo_symbol')]
  gene_mut_by_sample[, mut_type := fcase(has_snv & has_indel, 'SNV + Indel',
                                         has_snv, 'SNV',
                                         has_indel, 'Indel')]
  
  
  
  setnames(gene_mut_by_sample, 'hugo_symbol', 'gene')
  pathways_by_sample = merge.data.table(pway_by_gene, gene_mut_by_sample, by = 'gene',
                                        all = F, allow.cartesian = T)
  
  
  
  
  pw_mut_type = pathways_by_sample[path_id %in% pw_to_use, 
                                   .(mut_type = fcase(all(mut_type == 'SNV'), 'SBS', 
                                                      all(mut_type == 'Indel'), 'Indel',
                                                      default = 'SBS + Indel')), by = c('path_id', 'patient_id')]
  path_detail = dcast.data.table(pw_mut_type, patient_id ~ path_id, value.var = 'mut_type')
  

  filtered_genes_by_pathway = pway_by_gene[pw_to_use, on = 'path_id'][! duplicated(gene)]
  
  
  if(! is.null(gene_copy_calls)) {
    cna = fread(gene_copy_calls)[Unique_Patient_Identifier %in% included_samples]
    cna[filtered_genes_by_pathway, path_id := path_id, on = 'gene']
    path_cna = unique(cna[! is.na(path_id), .(path_id, patient_id = Unique_Patient_Identifier, copy)])
    
    cna_detail = path_cna[, .(cna_type = fcase(all(copy == 'neutral'), 'Neutral',
                                               all(copy %in% c('amp', 'neutral')), 'Copy gain',
                                               all(copy %in% c('del', 'neutral')), 'Copy decrease',
                                               default = 'Mixed')),
                          by = c('path_id', 'patient_id')]
    # We'll ignore the mixed CNA type
    combined_detail = merge.data.table(pw_mut_type, cna_detail, all = T, by = c('patient_id', 'path_id'))
  } else {
    combined_detail = copy(pw_mut_type)
    combined_detail[, cna_type := NA]
  }
  combined_detail[, merged_detail := fcase(cna_type %in% c('Neutral', 'Mixed') | is.na(cna_type), mut_type,
                                           cna_type == 'Copy gain' & ! is.na(mut_type), 'Copy gain + SBS/Indel',
                                           cna_type == 'Copy decrease' &  ! is.na(mut_type), 'Copy decrease + SBS/Indel',
                                           default = cna_type)]
  combined_detail[, patient_id := factor(patient_id, levels = included_samples)]
  output = dcast.data.table(combined_detail, patient_id ~ path_id, value.var = 'merged_detail', drop = FALSE)
  return(output)
}

# MSK samples are all IHC
msk_samples = cesa$samples[study == "MSK-IMPACT_iCCA", Unique_Patient_Identifier]
pw_in_msk = pathway_info[is_cancer_gene_path == TRUE, path_id]
msk_gene_copy_calls = 'gene_copy_calls/MSK_2021_IHCH_gene_copy.txt'

msk_landscape_features = prep_landscape_data(included_samples = msk_samples, pw_to_use = pw_in_msk, 
                                             gene_copy_calls = msk_gene_copy_calls)
fwrite(msk_landscape_features, 'output/landscape/MSK_path_cna_features.txt', sep = '\t', na = 'NA')

pan_landscape_features = prep_landscape_data(included_samples = cesa$samples$Unique_Patient_Identifier,
                                             pw_to_use = pathway_info$path_id)
fwrite(pan_landscape_features, 'output/landscape/pan_landscape_features.txt', sep = "\t", na = 'NA')

nontarget_samples = cesa$samples[coverage != 'targeted', Unique_Patient_Identifier]
exome_pan_landscape_features = prep_landscape_data(included_samples = nontarget_samples,
                                             pw_to_use = pathway_info$path_id)
fwrite(exome_pan_landscape_features, 'output/landscape/nontarget_pan_landscape_features.txt', sep = "\t", na = 'NA')

ihc_samples = cesa$samples[cca_type == 'IHC', Unique_Patient_Identifier]
ihc_landscape_features = prep_landscape_data(included_samples = ihc_samples,
                                             pw_to_use = pathway_info$path_id)
fwrite(ihc_landscape_features, 'output/landscape/IHC_landscape_features.txt', sep = "\t", na = 'NA')


ihc_nontarget = prep_landscape_data(included_samples = intersect(ihc_samples, nontarget_samples),
                                    pw_to_use = pathway_info$path_id)
fwrite(ihc_nontarget, 'output/landscape/nontarget_ihc_features.txt', sep = "\t", na = 'NA')


