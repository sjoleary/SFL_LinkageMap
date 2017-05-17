
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)


# load file
merge_stats <- readLines("results/LG12.merge")

# markers with constraints
map_extract <- as.data.frame(grep("Map", merge_stats, value = TRUE)) %>%
  rename(temp = `grep("Map", merge_stats, value = TRUE)`) %>%
  separate(temp, into = c("x", "y"), sep = "]") %>%
  separate(y, into = c("MAP", "LOCUS"), sep = ":") %>%
    select(MAP, LOCUS) %>%
    mutate(MAP = str_replace(MAP, "\"", " ")) %>%
    mutate(LOCUS = str_replace(LOCUS, "\"", " ")) %>%
    arrange(MAP)

# extract all RMSE values
map_rmse <- as.data.frame(grep("]" , merge_stats, invert = TRUE, value = TRUE)) %>%
    rename(temp = `grep("]", merge_stats, invert = TRUE, value = TRUE)`) %>%
    mutate(temp = str_trim(temp)) %>%
    droplevels()

# RMSE values for mean
rmse_mean <- map_rmse %>%
    filter(grepl("mean", temp)) %>%
    separate(temp, into = c("x", "CAT", "RMSE"), sep = " ") %>%
    mutate(RMSE = parse_number(RMSE)) %>%
    select(CAT, RMSE) %>%
    mutate(K = c(1:50))

ggplot(rmse_mean, aes(x = K, y = RMSE)) +
    geom_point() +
    theme_standard

# RMSE values for SD
rmse_sd <- map_rmse %>%
    filter(grepl("sd", temp)) %>%
    separate(temp, into = c("x", "y"), sep = 4) %>%
    separate(y, into = c("CAT", "RMSE"), sep = 3) %>%
    mutate(RMSE = parse_number(RMSE)) %>%
    mutate(CAT = str_trim(CAT)) %>%
    select(CAT, RMSE) %>%
    mutate(K = c(1:50))

ggplot(rmse_sd, aes(x = K, y = RMSE)) +
    geom_point() +
    theme_standard

# merge mean & sd
rmse <- left_join(rmse_mean, rmse_sd, by = "K") %>%
    rename(mean = RMSE.x,
           sd = RMSE.y) %>%
    select(K, mean, sd)

ggplot(rmse, aes(mean, sd)) +
     geom_point() + 
    theme_standard

# RMSE values map 1 & 2
rmse_map <- map_rmse %>%
    filter(!grepl("sd | mean", temp)) %>%
    filter(!grepl("map", temp)) %>%
    separate(temp, into = c("CAT", "y"), sep = c(1)) %>%
    separate(y, into = c("y", "RMSE"), sep = -6) %>%
    mutate(RMSE = parse_number(RMSE)) %>%
    mutate(CAT = str_trim(CAT)) %>%
    select(CAT, RMSE)

rmse_map1 <- rmse_map %>%
    filter(CAT == 1) %>%
    mutate(K = c(1:50))

ggplot(rmse_map1, aes(x = K, y = RMSE)) +
    geom_point() +
    theme_standard

rmse_map2 <- rmse_map %>%
    filter(CAT == 2) %>%
    mutate(K = c(1:50))

# create data frame with all rmse values
LG12_merge <- bind_rows(rmse_mean, rmse_sd)
LG12_merge <- bind_rows(LG12_merge, rmse_map1)
LG12_merge <- bind_rows(LG12_merge, rmse_map2)

ggplot(LG12_merge, aes(x = K, y = RMSE)) +
    geom_point() +
    facet_grid(CAT ~ ., scales = "free") +
    theme_facet

View(LG12[[22]])

    