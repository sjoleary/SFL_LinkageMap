
library(dplyr)
library(tidyr)

merge_stats <- readLines("results/LG10.merge")

map_extract <- as.data.frame(grep("Map", merge_stats, value = TRUE)) %>%
  rename(temp = `grep("Map", merge_stats, value = TRUE)`) %>%
  separate(temp, into = c("x", "y"), sep = "]") %>%
  separate(y, into = c("MAP", "LOCUS"), sep = ":")

map_rmse <- as.data.frame(grep("[1]", merge_stats, invert = TRUE, value = TRUE))
