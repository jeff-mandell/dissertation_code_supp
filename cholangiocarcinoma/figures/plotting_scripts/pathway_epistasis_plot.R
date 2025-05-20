library(data.table)
library(cancereffectsizeR)
library(stringr)
library(ggrepel)

# Plot pathway level epistasis effect
# Author: JL
# Last revised (WIP): 05/15/24, 01/16/25

comp_results = fread('./analysis/path_analysis/path_epistasis/epi_outputs_sig/path_epi_sig_out_panCCA.csv')

# Checking epistasis p-values:
comp_results$p_epistasis_original <- comp_results$p_epistasis
# All results are significant (even) with BH correction
comp_results$p_epistasis_bh <- p.adjust(comp_results$p_epistasis, method = 'BH') 
comp_results[p_epistasis_bh < 0.05, .N]
comp_results$p_epistasis <- comp_results$p_epistasis_original

comp_results <- comp_results[p_epistasis < 0.05,]

# comp_results$p_A_change_orig <- comp_results$p_A_change
# comp_results$p_A_change <- p.adjust(comp_results$p_A_change, method = 'BH')
# comp_results$p_B_change_orig <- comp_results$p_B_change
# comp_results$p_B_change <- p.adjust(comp_results$p_B_change, method = 'BH')

### DOT AND ARROW PLOT
## For filtering by mutational burden summed over all genes in pathway
# Same as in gene_epistasis_analysis.R:
cesa <- load_cesa("./output/cca_cesa.rds")

effectTb <- cesa$selection$all_effects
pway_by_variant = readRDS('./analysis/path_analysis/present_variants_by_path_cons_exome.rds')
# Function to sum included_with_variant for a given pathway
compute_sum <- function(variant_ids, effect_table) {
  sum(effect_table[variant_id %in% variant_ids, sum(included_with_variant, na.rm = TRUE)])
}

# Compute the sum for each pathway
pathway_mutation_counts <- data.table(
  PATH_ID = names(pway_by_variant),
  sum_included_with_variant = sapply(pway_by_variant, compute_sum, effect_table = effectTb)
)

# table(pathway_mutation_counts$sum_included_with_variant)
# 0    1    2    3    4    5    6    9 
# 3317   45    8    3    2    1    1    1 

# Filter to 16 pathways with at least 2 sum_included_with_variant
filtered_paths = pathway_mutation_counts[sum_included_with_variant >=2,]$PATH_ID
comp_results[pathway_mutation_counts, sum_included_with_variant.A := i.sum_included_with_variant, on=c('variant_A'='PATH_ID')]
comp_results[pathway_mutation_counts, sum_included_with_variant.A := i.sum_included_with_variant, on=c('variant_A'='PATH_ID')]

# Order by epistatic ratio
# Filter pathway pairs to those with significant path-level epistatic effect (i.e. p_A_change or p_B_change)
comp_results_plot <- comp_results[(p_A_change < .05) | (p_B_change < .05)] 

## For ordering by selection intensity (maximum effect reported from any cca type run)
cca_tb = fread('./analysis/path_analysis/pathway_effects.txt')

max_effectTb = as.data.table(unique(cca_tb$path_id))
colnames(max_effectTb) <- "path_id"
max_effectTb[cca_tb[cca_type == "IHC"], icca_selection_intensity := i.novel_variant_effect, on=c("path_id")]
max_effectTb[cca_tb[cca_type == "EHC"], ecca_selection_intensity := i.novel_variant_effect, on=c("path_id")]
max_effectTb[cca_tb[cca_type == "IHC"], dcca_selection_intensity := i.novel_variant_effect, on=c("path_id")]
max_effectTb[cca_tb[cca_type == "PHC"], pcca_selection_intensity := i.novel_variant_effect, on=c("path_id")]
max_effectTb[cca_tb[cca_type == "panCCA"], pancca_selection_intensity := i.novel_variant_effect, on=c("path_id")]
# some pathways don't have a corresponding row in a specific cca type run -> replace NA selection_intensities with 0.
max_effectTb[is.na(max_effectTb)] <- 0

# max selection intensities come from every cca type run, including pancca run
max_effectTb[, selection_intensity := max(dcca_selection_intensity, icca_selection_intensity, 
                                          pcca_selection_intensity, ecca_selection_intensity, 
                                          pancca_selection_intensity), by=c("path_id")]
# Get max and mean effect sizes for this output table
comp_results[max_effectTb, pathway_A_effect := i.selection_intensity, on=c("pathway_A_id" = "path_id")]
comp_results[max_effectTb, pathway_B_effect := i.selection_intensity, on=c("pathway_B_id" = "path_id")]
comp_results$max_effect <- apply(comp_results[, c("pathway_A_effect", "pathway_B_effect")], 1, max)
comp_results$mean_effect <- apply(comp_results[, c("pathway_A_effect", "pathway_B_effect")], 1, mean)

# Filter pathway pairs to those with significant path-level epistatic effect (i.e. p_A_change or p_B_change)
comp_results_plot <- comp_results[(p_A_change < .05) | (p_B_change < .05)] 

# Two methods to order the pathway pairs for plotting
# (A) Filter for top pathway pairs by mean effect size (of the 2 pathways)
comp_results_plot <- comp_results_plot[order(-mean_effect)]

# comp_results_plot <- unique( c(comp_results_plot$pathway_A_id,comp_results_plot$pathway_B_id) )
# 112 pairs

# unique_A_ids <- comp_results_plot[, .N, by = pathway_A_id][N == 1, pathway_A_id]
# unique_B_ids <- comp_results_plot[, .N, by = pathway_B_id][N == 1, pathway_B_id]
# comp_results_plot <- comp_results_plot[! (pathway_A_id %in% unique_B_ids),]
# comp_results_plot <- comp_results_plot[! (pathway_B_id %in% unique_A_ids),]

# comp_results_plot <- comp_results_plot[!duplicated(comp_results_plot$pathway_A_id), ] # 137 pairs
# comp_results_plot <- comp_results_plot[!duplicated(comp_results_plot$pathway_B_id), ] # 34 pairs
# comp_results_plot <- comp_results_plot[1:20, ]

# (B) Filter for top 10 pathway pairs by AB epistastic ratio
comp_results_plot <- comp_results_plot[order(-AB_epistatic_ratio)]

## Getting into plotting
# effectTb <- fread(file='./analysis/path_analysis/path_effects/path_panCCA_effectTb.txt', sep='\t')
# comp_results_plot[effectTb, ces_A_null := i.selection_intensity, on=c(pathway_A_id = "pathway_id")]
# comp_results_plot[effectTb, ces_B_null := i.selection_intensity, on=c(pathway_B_id = "pathway_id")]

## Pretty print pathway names -- WIP --
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("reactome.db")

# From Bioconductor reactome.db package:
# xx <- as.data.table(reactome.db::reactomePATHID2NAME)
# xx <- xx[grepl("^Homo sapiens:*", xx$path_name), ]
# xx$path_name <- gsub("^Homo sapiens: *", "", xx$path_name)
# xx[, path_name_upper := toupper(path_name)]

path_labelmaker <- function(x) {
  path_name_lst <- unlist(strsplit(x, ", "))
  path_label = list()
  for (pth in path_name_lst){
    path_name <- str_remove_all(gsub("\\..*","",pth),"_PATHWAY")
    if (str_detect(path_name, "WP")){
      path_name <- stringr::str_to_sentence(substr(path_name, 4, nchar(path_name)))
    }
    else {
      path_name <- stringr::str_to_sentence(path_name)
    }
    path_label <- append(path_label, path_name)
  }
  path_label_final <- paste0(path_label, collapse=", ")
  return(path_label_final)
}

comp_results_plot[, original_A := variant_A]
comp_results_plot[, original_B := variant_B]

comp_results_plot[, plot_name_long_A := sapply(comp_results_plot$original_A, path_labelmaker)]
comp_results_plot$plot_name_long_A <- gsub("_", " ", comp_results_plot$plot_name_long_A)
comp_results_plot[, variant_A := sapply(plot_name_long_A, function(x) paste0(strwrap(x, width=40), collapse="\n"))]

comp_results_plot[, plot_name_long_B := sapply(comp_results_plot$original_B, path_labelmaker)]
comp_results_plot$plot_name_long_B <- gsub("_", " ", comp_results_plot$plot_name_long_B)
comp_results_plot[, variant_B := sapply(plot_name_long_B, function(x) paste0(strwrap(x, width=40), collapse="\n"))]

# comp_results_plot[xx, variant_A := i.path_name, on=c("original_A" = "path_name_upper")]
# comp_results_plot[xx, variant_B := i.path_name, on=c("original_B" = "path_name_upper")]

# Plot desired pathway/variant pairs
plot_epistasis(comp_results_plot, pairs_per_row = 6, variant_label_size = 6.5)

ggsave(file='./figures/pathway_plots/path_epistasis_eCCA_top10MeanEffect.png')
