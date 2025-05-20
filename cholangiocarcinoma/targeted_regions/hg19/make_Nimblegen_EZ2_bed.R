library(data.table)

intervals = fread('NimbleGenEZ2Tiled_hg19_44Mb_intervals.txt', header = F)
intervals[, c("chr", "start", "end") := tstrsplit(V1, split = '[:-]')]
gr = makeGRangesFromDataFrame(intervals)
gr = reduce(sort(gr))
export(gr, 'NimbleGen_EZ2_hg19.bed')
