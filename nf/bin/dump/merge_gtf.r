args <- commandArgs(trailingOnly = TRUE)

library(rtracklayer)
library(plyranges)

human_gtf_path <- args[1]
sinv_gtf_path <- args[2]
output_path <- args[3]

setwd(dirname(output_path))

human_GTF <- import(human_gtf_path)
SINV_GTF <- import(sinv_gtf_path)

seqlevels(SINV_GTF) <- "SINV"
SINV_GTF <- SINV_GTF %>% mutate(seqnames = "SINV")

export(SINV_GTF, "test.gtf", "gtf")

human_SINV_GTF <- GRangesList(human_GTF, SINV_GTF) %>% unlist()

#export(human_SINV_GTF, output_path)

