library(tidyr)
library(dplyr)
library(stringr)

SNPs_filtered <- read.table("data/VCF/FamB.SNP.ldepth.mean",
                            header = TRUE, stringsAsFactors = FALSE) %>%
    select(CHROM, POS)

mapped_markers <- read.table("results/FamA_family.map", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(marker = (str_sub(LOCUS, 2))) %>%
    select(marker)

mapped_SNPs <- filter(SNPs_filtered, CHROM %in%mapped_markers$marker)

