##### FUNCTIONS TO FORMAT & PLOT SYNTENY MAPPER RESULTS FILES #####

## function to read in *synteny.stats files ----------------------- 

# load libraries 
library(dplyr) 
library(plyr) 

# arguments: 
# dir (directory where stats files are located) 
# species (vector with species abbreviations) 

# columns: 
# N_LOCI = Total Loci in Blocks 
# N_BLOCKS = Total Blocks 
# MEAN_LOCI = Mean loci per block 
# MIN_LOCI = Min loci per block 
# MAX_LOCI = Max loci per block 
# TOTAL_BP = Total Comp Block Size (bp) 
# MEAN_BP = Mean Comp Block Size (bp) 
# MIN_BP = Min Comp Block Size (bp) 
# MAX_BP = Max Comp Block Size (bp) 
# TOTAL_cM = Total Map Block Size (cM) 
# MEAN_cM = Mean Map Block Size (cM) 
# MIN_cM = Min Map Block Size (cM) 
# MAX_cM = Max Map Block Size (cM) 

read.synteny.stats <- function(dir, species) { 
  temp <- lapply(species, function(sp){ 
    filename <- paste(sp, ".synteny.stats", sep = "") 
    path <- file.path(dir, filename) 
    synteny.stats <- read.table(path, header = FALSE, sep = ":", 
                                col.names = c("STAT", "species")) 
    synteny.stats <- setNames(data.frame(t(synteny.stats[, -1])), synteny.stats[,1]) %>% 
      mutate(SPECIES = sp) %>% 
      select(SPECIES, 1:13) 
  }) 
  ldply(temp, rbind) %>% 
    dplyr::rename(N_LOCI = `Total Loci in Blocks`, 
           N_BLOCKS = `Total Blocks`, 
           MEAN_LOCI = `Mean loci per block`, 
           MIN_LOCI = `Min loci per block`, 
           MAX_LOCI = `Max loci per block`, 
           TOTAL_BP = `Total Comp Block Size (bp)`, 
           MEAN_BP = `Mean Comp Block Size (bp)`, 
           MIN_BP = `Min Comp Block Size (bp)`, 
           MAX_BP = `Max Comp Block Size (bp)`, 
           TOTAL_cM = `Total Map Block Size (cM)`, 
           MEAN_cM = `Mean Map Block Size (cM)`, 
           MIN_cM = `Min Map Block Size (cM)`, 
           MAX_cM = `Max Map Block Size (cM)`) 
}   

## function to plot synteny stats overview ------------------------ 
# expects synteny.stats file in format of read.synteny.stats

# libraries
library(tidyr)
library(ggplot2)

plot.synteny.stats <- function(synteny_stats) { 
  tidy <- synteny_stats %>% 
    gather("STAT", "RESULT", 2:14) 
  ggplot(tidy, aes(x = SPECIES, y = RESULT, fill = SPECIES)) + 
    geom_bar(color = "black", stat = "identity") + 
    facet_wrap( ~ STAT, scales = "free") 
  
} 

## function to plot loci per block overview ----------------------- 
# expects synteny.stats file in format of read.synteny.stats

# libraries
library(tidyr)
library(ggplot2)
library(ggthemes)

plot.synteny.loci <- function(synteny_stats) { 
  
  # Format data set
  tidy <- synteny_stats %>% 
    gather("STAT", "RESULT", 2:14) %>%
    filter(STAT %in% c("MEAN_LOCI", "MIN_LOCI", "MAX_LOCI"))
  
  # plot
  ggplot(tidy, aes(x = SPECIES, y = RESULT, fill = SPECIES)) + 
    geom_bar(data = filter(tidy, STAT == "MEAN_LOCI"),
             color = "black", stat = "identity") +
    geom_point(data = filter(tidy, STAT == "MIN_LOCI"),
               shape = 21, size = 4, fill = "black") +
    geom_point(data = filter(tidy, STAT == "MAX_LOCI"),
               shape = 21, size = 4, fill = "black") +
    labs(x = "Species genome", y = "Mean number of loci (Max/Min)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      strip.background = element_rect(fill = "grey95", color = "black"), 
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
} 

## function to plot total number of loci, blocks, length ---------- 
# expects synteny.stats file in format of read.synteny.stats

# libraries
library(tidyr)
library(ggplot2)

plot.synteny.total <- function(synteny_stats) { 
  
  # Format data set
  tidy <- synteny_stats %>% 
    gather("STAT", "RESULT", 2:14) %>%
    filter(STAT %in% c("N_LOCI", "N_BLOCKS", "TOTAL_BP", "TOTAL_cM"))
  
  # plot
  ggplot(tidy, aes(x = SPECIES, y = RESULT, fill = SPECIES)) + 
    geom_bar(color = "black", stat = "identity") +
    labs(x = "Species genome", y = "Total number/length") +
    facet_wrap(~ STAT, scales = "free") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA), 
      panel.border = element_rect(fill = NA, color = "black"), 
      panel.grid.major = element_line(color = "grey85"), 
      panel.grid.minor = element_blank(), 
      strip.background = element_rect(fill = "grey95", color = "black"), 
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))  
} 

## function to plot length of syntenic blocks on chromosome (bp) ----
# expects synteny.stats file in format of read.synteny.stats

# libraries
library(tidyr)
library(ggplot2)

# plot bp overview
plot.synteny.bp <- function(synteny_stats) { 
  
  # Format data set
  tidy <- synteny_stats %>% 
    gather("STAT", "RESULT", 2:14) %>%
    filter(STAT %in% c("MEAN_BP", "MIN_BP", "MAX_BP"))
  
  # plot
  ggplot(tidy, aes(x = SPECIES, y = RESULT, fill = SPECIES)) + 
    geom_bar(data = filter(tidy, STAT == "MEAN_BP"),
             color = "black", stat = "identity") +
    geom_point(data = filter(tidy, STAT == "MIN_BP"),
               shape = 21, size = 4, fill = "black") +
    geom_point(data = filter(tidy, STAT == "MAX_BP"),
               shape = 21, size = 4, fill = "black") +
    labs(x = "Species genome", y = "Mean length of synteny blocks in bp (Max/Min)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      strip.background = element_rect(fill = "grey95", color = "black"), 
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
} 

## function to plot length of syntenic blocks on LG (cM) ---------- 
# expects synteny.stats file in format of read.synteny.stats

# libraries
library(tidyr)
library(ggplot2)

plot.synteny.cM <- function(synteny_stats) { 
  
  # Format data set
  tidy <- synteny_stats %>% 
    gather("STAT", "RESULT", 2:14) %>%
    filter(STAT %in% c("MEAN_cM", "MIN_cM", "MAX_cM"))
  
  # plot
  ggplot(tidy, aes(x = SPECIES, y = RESULT, fill = SPECIES)) + 
    geom_bar(data = filter(tidy, STAT == "MEAN_cM"),
             color = "black", stat = "identity") +
    geom_point(data = filter(tidy, STAT == "MIN_cM"),
               shape = 21, size = 4, fill = "black") +
    geom_point(data = filter(tidy, STAT == "MAX_cM"),
               shape = 21, size = 4, fill = "black") +
    labs(x = "Species genome", y = "Mean length of synteny blocks in bp (Max/Min)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      strip.background = element_rect(fill = "grey95", color = "black"), 
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
} 

## function to read in all.blocks.tab -----------------------------
# argument: dir = path to data folder

read.all.blocks <- function(dir) { 
  path <- file.path(dir, "all.blocks.tab") 
  synt_blocks <- read.table(path, header = TRUE,
                            colClasses = c("factor", "factor", "factor", "numeric", 
                                           "factor", "integer", "integer", "integer",
                                           "factor", "numeric", "numeric", "numeric"))
  # order MAP_LG factor
  synt_blocks$MAP_LG <- factor(synt_blocks$MAP_LG, 
                               levels = c("1", "2", "3", "4", "5", "6", "7", "8", 
                                          "9", "10", "11", "12", "13", "14", "15", "16",
                                          "17", "18", "19", "20", "21", "22", "23", "24"))
  synt_blocks
}   

## function to create Oxford bubble plots --------------------------
# expects synteny.blocks file in format of read.all.blocks
# arguments 
  # synt_blocks = data.frame in format all.blocks.tab
  # species = vector of species to be plotted

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)

plot.oxford <- function(synt_blocks, species) {
  # identify homologous linkage groups
  Homolog <- synt_blocks %>%
    group_by(COMP_SPECIES) %>%
    arrange(MAP_LG) %>%
    select(COMP_SPECIES, MAP_LG, COMP_CHR) %>%
    unite(HOMOLOG, 1:3, sep = ">") %>%
    distinct() %>%
    separate(HOMOLOG, into = c("COMP_SPECIES", "MAP_LG", "COMP_CHR"), sep = ">")
  
  # order Linkage groups
  Homolog$MAP_LG <- factor(Homolog$MAP_LG, 
                           levels = c("1", "2", "3", "4", "5", "6", "7", "8", 
                                      "9", "10", "11", "12", "13", "14", "15", "16",
                                      "17", "18", "19", "20", "21", "22", "23", "24"))
  
  lapply(species, function(species){
    
    # order Chromosomes according to LG
    homologs <- Homolog %>%
      filter(COMP_SPECIES %in% species) %>%
      distinct(MAP_LG) %>%
      arrange(MAP_LG)
    CHR <- as.factor(homologs$COMP_CHR)
    
    # determine number of blocks per chromosome/linkage group combinations
    temp <- synt_blocks %>%
      filter(COMP_SPECIES %in% species) %>%
      group_by(MAP_LG) %>%
      count(MAP_LG, COMP_CHR) %>%
      filter(!is.na(COMP_CHR))
    temp$COMP_CHR <- factor(temp$COMP_CHR, levels = unique(CHR))
    # plot LG vs CHR and number of blocks
    ggplot(temp, aes(x = MAP_LG, y = COMP_CHR)) +
      geom_point(aes(size = n), 
                 shape = 21, color = "black", fill = "grey85") +
      scale_size_continuous(range=c(2,15)) +
      labs(x = "Linkage Group Southern Flounder", y = paste("Chromosome", species)) +
      theme_standard
  })
}

## function to dataframe of corresponding LG/CHR --------------------
# expects synteny.blocks file in format of read.all.blocks

# libraries
library(dplyr)

LG.Homologs <- function(synteny_blocks){
  synt_blocks %>%
    group_by(COMP_SPECIES) %>%
    arrange(MAP_LG) %>%
    select(COMP_SPECIES, MAP_LG, COMP_CHR) %>%
    unite(HOMOLOG, 1:3, sep = ">") %>%
    distinct() %>%
    separate(HOMOLOG, into = c("COMP_SPECIES", "MAP_LG", "COMP_CHR"), sep = ">")
}