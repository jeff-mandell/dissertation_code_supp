library(data.table)

# All patient samples are in the sample key. And sample are used in the CESAnalysis.
sk = fread('combined_sample_key.txt')
sk[Unique_Patient_Identifier %like% 'nakamura', study := 'nakamura']
sk[Unique_Patient_Identifier %like% 'jusakul', 
   study := paste0('jusakul.', coverage)]

# One sequencing type per study, with the above tweak
stopifnot(uniqueN(sk[, .(study, coverage)]) == uniqueN(sk$study))
counts = dcast(sk[, .N, by = c('study', 'cca_type')], 
               study ~ cca_type, fill = 0, value.var = 'N')
counts[, total := IHC + PHC + DCC + EHC]
setcolorder(counts, c('study', 'IHC', 'PHC', 'DCC', 'EHC', 'total'))

counts = counts[, s2 := tolower(study)][order(s2)][, s2 := NULL][]
counts = rbind(counts[study == 'yale'], counts[study != 'yale'])

fwrite(counts, 'output/study_included_sample_counts.txt', sep = "\t")

counts[, lapply(.SD, sum), .SDcols =  is.numeric]
# IHC   PHC   DCC   EHC total
# <int> <int> <int> <int> <int>
#   1:   930   186    74    15  1205



