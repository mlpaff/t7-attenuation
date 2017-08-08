library(tidyverse)
library(cowplot)

rna <- read.csv("rna_abundance.csv")

# wt background - genes 9 and 10, separated by replicate
rna %>% group_by(strain, rep) %>% 
  filter(strain %in% c('T7Hi', 'phi9v2', 'phi10', 'phi910v2'), 
         grepl('p42|p43', gene)) %>% 
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', "gene 10")) %>%
  ggplot(aes(x=strain, y=rel_ab, fill=strain)) + geom_bar(position='dodge', stat='identity') + facet_wrap(~gene + rep)

# plot mean rel. abundance for genes 9 and 10 
rna %>% group_by(strain) %>% 
  filter(strain %in% c('T7Hi', 'phi9v2', 'phi10', 'phi910v2'), 
         grepl('p42|p43', gene)) %>%
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', 'gene 10')) %>% group_by(strain, gene) %>% 
  summarize(mean=mean(rel_ab), se=sd(rel_ab)/sqrt(n())) %>%
  ggplot(aes(x=strain, y=mean, fill=strain)) + geom_bar(position='dodge', stat='identity') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4, position='dodge', stat='identity') + facet_wrap(~gene)

rna %>% group_by(strain, rep) %>%
  filter(strain %in% c('44', '44-phi9v2', '44-phi10'), 
         grepl('p42|p43', gene)) %>% 
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', "gene 10")) %>%
  ggplot(aes(x=strain, y=rel_ab, fill=strain)) + geom_bar(position='dodge', stat='identity') + facet_wrap(~gene + rep)

rna %>% group_by(strain) %>%
  filter(strain %in% c('44', '44-phi9v2', '44-phi10'), 
         grepl('p42|p43', gene)) %>%
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', 'gene 10')) %>% group_by(strain, gene) %>% 
  summarize(mean=mean(rel_ab), se=sd(rel_ab)/sqrt(n())) %>%
  ggplot(aes(x=strain, y=mean, fill=strain)) + geom_bar(position='dodge', stat='identity') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4, position='dodge', stat='identity') + facet_wrap(~gene)

rna %>% group_by(strain, rep) %>%
  filter(strain %in% c('42', '42-phi9v2', '42-phi10'), 
         grepl('p42|p43', gene)) %>% 
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', "gene 10")) %>%
  ggplot(aes(x=strain, y=rel_ab, fill=strain)) + geom_bar(position='dodge', stat='identity') + facet_wrap(~gene + rep)

rna %>% group_by(strain, rep) %>%
  filter(strain %in% c('T7Hi', 'phi9st', 'phi910st'), 
         grepl('p42|p43', gene)) %>% 
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', "gene 10")) %>%
  ggplot(aes(x=strain, y=rel_ab, fill=strain)) + geom_bar(position='dodge', stat='identity') + facet_wrap(~gene + rep)

rna %>% group_by(strain) %>%
  filter(strain %in% c('T7Hi', 'phi9st', 'phi910st'),  
         grepl('p42|p43', gene)) %>%
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', 'gene 10')) %>% group_by(strain, gene) %>% 
  summarize(mean=mean(rel_ab), se=sd(rel_ab)/sqrt(n())) %>%
  ggplot(aes(x=strain, y=mean, fill=strain)) + geom_bar(position='dodge', stat='identity') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4, position='dodge', stat='identity') + facet_wrap(~gene)

# All combined
rna %>% group_by(strain, rep) %>% 
  filter(grepl('p42|p43', gene)) %>%
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', "gene 10")) %>%
  ggplot(aes(x=strain, y=rel_ab, fill=strain)) + geom_bar(position='dodge', stat='identity') + facet_wrap(~gene + rep)

rna %>% group_by(strain) %>%
  filter(grepl('p42|p43', gene)) %>%
  mutate(rel_ab=tpm/1000000, gene=if_else(length==1611, 'gene 9', 'gene 10')) %>% group_by(strain, gene) %>% 
  summarize(mean=mean(rel_ab), se=sd(rel_ab)/sqrt(n())) %>%
  ggplot(aes(x=strain, y=mean, fill=strain)) + geom_bar(position='dodge', stat='identity') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4, position='dodge', stat='identity') + facet_wrap(~gene)
