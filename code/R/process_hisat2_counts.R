# Combine read-counts data generated from hisat2 and normalize data
require(readr)
require(purrr)
library(stringr)
library(openxlsx)
library(tidyverse)
library(cowplot)

# Path to directory containing counts files 
dir_path <- "../../data/hisat2_data/"

# find and list path to all '*counts.txt' files
files <- dir(dir_path, pattern='*.txt')
# list sample numbers
sample_list <- str_match(files, '(PA[0-9]+)')[,1]

# Information on strain, treatment, genetic background, and promoter knockout type to be added to the data frame of raw RNA counts
conditions <-data.frame(
  strain=c('T7Hi', 'phi9v2', 'phi910v2','phi10', '11-44', '11-44-phi9v2', '11-44-phi10', '8st-9', '8st-910', '910L2evo', '8st-910evo'),
  treatment=c('wt', 'wt-9', 'wt-910', 'wt-10', 'deop', 'deop-9', 'deop-10', '8st-9', '8st-910', '910L2evo', '8st-910evo'),
  background=c('wt', 'wt', 'wt', 'wt', '10deop', '10deop', '10deop', '8st', '8st', 'wt', '8st'),
  knockout=c('wt', 'phi9', 'phi910', 'phi10', 'wt', 'phi9', 'phi10', 'phi9', 'phi910', 'phi910', 'phi910')
)

# read in sample list
samples <- read.xlsx("../../data/sample_list.xlsx") %>% select(Sample, strain=Strain, rep=Biological.replicate) %>%
  # Filter samples from the sample_list
  filter(Sample %in% sample_list) %>% 
  left_join(data_frame(path=files, Sample=sample_list)) 

colnames <- c("ref_genome", "start", "stop", "locus_tag", "gene", "counts")

# read in counts data
data <- samples %>%
  mutate(file_contents=map(path, ~ read_table2(file.path(dir_path, .), col_names=colnames, col_types=cols()))) %>%
  unnest() %>% select(-path)

rna <- left_join(conditions, data)

# Normalize the raw RNA counts, first to the sum of raw counts for all genes up to and including gene 7.7, then by TPM
rna %>% group_by(rep, strain) %>% 
  mutate(nfactor=sum(counts[stop <= 20240])) %>%
  mutate(normcounts=counts/nfactor) %>%
  mutate(rpk=normcounts/((stop-start)/1000)) %>%
  mutate(rpm=sum(rpk)) %>%
  mutate(tpm=rpk/rpm) -> counts

counts$gene <- factor(counts$gene, levels=rna$gene[1:60])

write.csv(counts, "../../data/results/counts_rna_abundance.csv", row.names=FALSE)

rna %>% filter(rep == 1, strain == 'T7Hi', stop <= 20240)
