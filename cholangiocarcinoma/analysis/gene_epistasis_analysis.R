# Purpose: Run pairwise epistasis for each unique gene from panCCA selection output.

library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg38)

# variants and their unique genes
# (1) Consensus exome run
# (2) Total consensus run
# Compare 
cesa <- load_cesa("./output/cca_cesa.rds")

# For each pair, if BOTH genes are covered in targeted MSK (see list), use  exome x MSK variants (i.e. variants covered by both) to get genes to use, which allows us to increase the samples used in the inference. 
# 307 genes covered by both exome and MSK.
msk_inter = Reduce(intersect, cesa$coverage_ranges$targeted[c("MSK341_ip100", "MSK410_ip100", "MSK468_ip100")])
exome_inter = Reduce(intersect, cesa$coverage_ranges$exome)
msk_exome_inter = intersect(msk_inter, exome_inter)

variants_in_inter =  select_variants(cesa, gr = msk_exome_inter)[(variant_type == 'snv' & essential_splice == TRUE) | 
                                                                   (variant_type == 'aac' & aa_ref != aa_alt)]
msk_exome_genes = variants_in_inter[, .(N = sum(maf_prevalence)), by = 'gene'][N > 19, gene] # 29 genes
msk_exome_vars_for_epistasis = variants_in_inter[gene %in% msk_exome_genes]

# 75 genes covered only in exome, but still run msk_exome_genes in exome run (just won't use all the samples) for when only 1 of the genes is in the msk_exome_genes list
# Note: Later, we can filter out the duplicate runs
exome_variants =  select_variants(cesa, gr = exome_inter)[(variant_type == 'snv' & essential_splice == TRUE) | 
                                                                   (variant_type == 'aac' & aa_ref != aa_alt)]
exome_genes = exome_variants[, .(N = sum(maf_prevalence)), by = 'gene'][N > 19, gene] # 70T genes
# all(msk_exome_genes %in% exome_genes) # TRUE
exome_vars_for_epistasis = exome_variants[gene %in% exome_genes]

# Run pan and subtypes independently
# Run panCCA epistasis
print("Running epistasis for panCCA")
# 29 genes in 10 mins
cesa = ces_gene_epistasis(cesa = cesa, 
                          genes = msk_exome_genes, 
                          variants = msk_exome_vars_for_epistasis, 
                          samples = cesa$samples, 
                          run_name = 'msk_exome_panCCA')
msk_exome_output_panCCA = cesa$epistasis$msk_exome_panCCA
# 70 genes in 1hr 10mins
cesa = ces_gene_epistasis(cesa = cesa, 
                          genes = exome_genes, 
                          variants = exome_vars_for_epistasis, 
                          samples = cesa$samples, 
                          run_name = 'exome_panCCA')
exome_output_panCCA = cesa$epistasis$exome_panCCA
fwrite(msk_exome_output_panCCA, file='./output/epistasis_output/gene_epi_MSKexome_panCCA.csv')
fwrite(exome_output_panCCA, file='./output/epistasis_output/gene_epi_exome_panCCA.csv')

# Don't run subtype epistasis yet...
# iCCA, dCCA, pCCA
cca_types = c('IHC', 'DCC', 'PHC')
for (this_cca_type in cca_types){
  print(paste('Running epistasis for:', this_cca_type))
  msk_run_name = paste0('msk_exome_', this_cca_type)
  exome_run_name = paste0('exome_', this_cca_type)
  cesa = ces_gene_epistasis(cesa = cesa, 
                            genes = msk_exome_genes, 
                            variants = msk_exome_vars_for_epistasis, 
                            samples = cesa$samples[cca_type == this_cca_type], 
                            run_name = msk_run_name)
  msk_exome_output = cesa$epistasis[[msk_run_name]]
  cesa = ces_gene_epistasis(cesa = cesa, 
                            genes = exome_genes, 
                            variants = exome_vars_for_epistasis, 
                            samples = cesa$samples[cca_type == this_cca_type], 
                            run_name = exome_run_name)
  exome_output = cesa$epistasis[[exome_run_name]]
  msk_filename = paste0('./output/epistasis_output/gene_epi_MSKexome_', this_cca_type, '.csv')
  exome_filename = paste0('./output/epistasis_output/gene_epi_exome_', this_cca_type, '.csv')
  fwrite(msk_exome_output, file=msk_filename)
  fwrite(exome_output, file=exome_filename)
}

# eCCA
print("Running epistasis for eCCA")
cesa = ces_gene_epistasis(cesa = cesa, 
                          genes = msk_exome_genes, 
                          variants = msk_exome_vars_for_epistasis, 
                          samples = cesa$samples[cca_type %in% c('EHC', 'DCC', 'PHC')], 
                          run_name = 'msk_exome_EHC')
cesa = ces_gene_epistasis(cesa = cesa, 
                          genes = exome_genes, 
                          variants = exome_vars_for_epistasis, 
                          samples = cesa$samples[cca_type %in% c('EHC', 'DCC', 'PHC')], 
                          run_name = 'exome_EHC')
msk_exome_output_eCCA = cesa$epistasis$msk_exome_EHC
exome_output_eCCA = cesa$epistasis$exome_EHC
fwrite(msk_exome_output_eCCA, file='./output/epistasis_output/gene_epi_MSKexome_EHC.csv')
fwrite(exome_output_eCCA, file='./output/epistasis_output/gene_epi_exome_EHC.csv')
