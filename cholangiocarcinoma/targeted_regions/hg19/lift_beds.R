library(rtracklayer)


beds = list.files(pattern = 'bed$')

chain_file = '~/reference/chains/hg19ToHg38.over.chain'
chain = import.chain(chain_file)
names(chain) = sub("^chr", "", names(chain))

for(bed in beds) {
  gr = import.bed(bed)
  seqlevelsStyle(gr) = 'NCBI'
  lifted = sort(reduce(unstrand(unlist(liftOver(gr, chain)))))
  seqlevelsStyle(lifted) = 'NCBI'
  
  prop = sum(width(lifted)) / sum(width(reduce(gr)))
  if(prop > 1 || prop < .95) {
    warning(bed, ': ', prop)
  }
  output = sub('\\.bed', '_hg38.bed.gz', bed)
  export.bed(lifted, output)
}
