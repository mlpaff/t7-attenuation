# Process post-Kallisto RNA-seq alignment data
require(tidyr)
require(readr)
require(purrr)
require(dplyr)
library(openxlsx)

# Path to directory containing RNA abundance tables 
dir_path <- "../../data/kallisto_output/"
# find and list path to all 'abundance.tsv' files
files <- dir(dir_path, pattern='*.tsv', recursive=TRUE)
# list of all *.fastq.gz files -- specifies which samples will be filtered from the master sample list
fastq_list <- paste0(basename(list.dirs(dir_path, full.names = FALSE, recursive = FALSE)), '.fastq.gz')

# read in sample list
samples <- read.xlsx("../../data/sample_list.xlsx") %>% select(Sample, Experiment.group, Strain, Biological.replicate, File) %>%
  rename(strain="Strain", rep="Biological.replicate") %>%
  # Filter samples for which there is a corresponding abundance.tsv file
  filter(File %in% fastq_list) %>% 
  left_join(data_frame(path=files, File=fastq_list))

# read in anundance data
data <- samples %>%
  mutate(file_contents=map(path, ~ read_tsv(file.path(dir_path, .), col_types=cols()))) %>%
  unnest() %>% select(-path, -File, -Experiment.group, -Sample) %>%
  # separate information from "target_id" into individual columns
  separate(target_id, c("organism", "gene", "prot_id"), sep='\\|', extra='drop')

write.csv(data, "../../data/results/kallisto_rna_abundance.csv")
