setwd("/home2/2459810a/iCLIP_analysis/NOP2_Innes/")
library(tidyverse)

read_coverage <- function(x){
  name <- basename(x)
  read_tsv(x, col_names = FALSE, col_types = "cnnn") %>%
    mutate(replicate = name)
}

get_coverage <- function(sample){
  IP_files <- list.files(dir, pattern = paste( "NOP2_IP_", sample, "*", sep = ""), full.names = TRUE)
  IP_files <- IP_files[c(2, 4, 6)]
  print(IP_files)
  IP_coverage <- map_df(IP_files, read_coverage) %>%
    filter(X1 == "SINV") %>%
    mutate(condition = "IP")
  SMI_files <- list.files(dir, pattern = paste( "NOP2_SMI_",sample, "*", sep = ""), full.names = TRUE)
  SMI_files <- SMI_files[c(2, 4, 6)]
  print(SMI_files)
  SMI_coverage <- map_df(SMI_files, read_coverage) %>%
    filter(X1 == "SINV") %>%
    mutate(condition = "SMI")
  df <- bind_rows(IP_coverage, SMI_coverage) %>% add_column(sample = sample)
  return(df)
}


##
dir <- "./Xlsite/shifted/bedgraph"

mock <- get_coverage("mock")
SINV <- get_coverage("SINV")


df <- rbind(mock, SINV)

write.csv(df, "GENOMEDIR/NOP2_SINV_coverage.csv")
