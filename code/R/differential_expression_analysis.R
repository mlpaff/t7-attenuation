# Differential expression analysis 
library(tidyverse)
library(cowplot)

tpm_stats <- read.csv("../../data/results/pairwise_t_vs_wt.csv")
deseq2_stats <- read.csv("../../data/results/deseq2_results_vs_wt.csv")

tpm_stats %>% filter(knockout=="phi9") #, gene %in% c('8', '9', '10A', '11', '12', '13'))
tpm_stats %>% filter(strain=='11-44', control=='T7Hi')
tpm_stats %>% filter(control=="11-44")
tpm_stats %>% filter(knockout=='phi10', control=="T7Hi")
tpm_stats %>% filter(knockout=='phi910', control=="T7Hi", gene %in% c('8', '9', '10A', '11', '12', '13'))

tpm_stats %>% filter(knockout=='phi910', control %in% c("8st-910", "phi910v2"), gene %in% c('8', '9', '10A', '11', '12', '13'))

deseq2_stats %>% filter(knockout=="phi9", gene %in% c('8', '9', '10A', '11', '12', '13'))
