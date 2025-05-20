library(cancereffectsizeR)

# The targeted regions files for most exome capture kits used in the analysis use hg19 coordinates.
# This script reads in hg19 BED files, coverts coordinates to hg38, and prints new files.

beds = list.files(path = 'targeted_regions/hg19/', pattern = 'bed.gz$', full.names = TRUE)
chain_file = 'reference/chains/hg19ToHg38.over.chain'
chain = import.chain(chain_file)

for (bed in beds) {
  outfile = paste0('targeted_regions/', basename(bed))
  outfile = sub('(_hg19)?\\.bed\\.gz', '_hg38.bed.gz', outfile)
  lift_bed(bed = bed, chain = chain, outfile = outfile)
}

