# Analyzing RNA count data from Hisat2 alignments, normalizing by calculating TPM 
library(tidyverse)
library(cowplot)

cbPalette <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999')

rna <- read.csv("../../data/restults/hisat2_rna_abundance.csv")

rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)
rna$strain <- factor(rna$strain, 
                     levels=c('T7Hi', 'phi9v2', 'phi10', 'phi910v2', '910L2evo', 
                              '11-44', '11-44-phi9v2', '11-44-phi10', 
                              '8st-9', '8st-910', '8st-910evo'))
rna$knockout <- factor(rna$knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'), 
                       labels=c(
                         'wt' = 'wt',
                         'phi9' = expression(paste(Delta, phi, '9')),
                         'phi10' = expression(paste(Delta, phi, '10')),
                         'phi910' = expression(paste(Delta, phi, '9/10'))))
rna$background <- factor(rna$background, 
                         levels=c('wt', '10deop', '8st'))
strain_labs <- c(
  'T7Hi' = 'wt',
  'phi9v2' = expression(paste(Delta, phi, '9'['wt'])),
  'phi910v2' = expression(paste(Delta, phi, '9/10'['wt'])),
  'phi10' = expression(paste(Delta, phi, '10'['wt'])),
  '11-44' = expression(paste(Delta, phi, '10'['deop'])),
  '11-44-phi9v2' = expression(paste(Delta, phi, '9'['10'['deop']])),
  '11-44-phi10' = expression(paste(Delta, phi, '10'['10'['deop']])),
  '8st-9' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])),
  '8st-910' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']])),
  '910L2evo' = expression(paste(Delta, phi, '9/10evo'['wt-L2'])),
  '8st-910evo' = expression(paste(Delta, phi, '9/10evo'['8'[Delta]['stop']]))
)
back_labs <- c(
  'wt' = 'wt',
  '8st'=expression(paste('8'[Delta]['stop'])),
  '10deop'=expression(paste('10'['deop']))
)

rna %>% group_by(rep, strain) %>% mutate(rpk=counts/((stop-start)/1000)) %>% 
  mutate(tpm=rpk/sum(rpk)/1000000) -> counts

counts$gene <- factor(counts$gene, levels=rna$gene[1:60], labels=paste('gene', rna$gene[1:60]))

# Plotting counts for genes 9 and 10A in the wt background. 
counts %>% filter(strain %in% c('T7Hi', 'phi9v2', 'phi10', 'phi910v2'), gene %in% c('gene 9', 'gene 10A'), rep!='2') %>% 
  group_by(strain) %>%
  ggplot(aes(x=strain, y=tpm, fill=strain)) + 
  stat_summary(geom="bar", fun.data='mean_se') + 
  geom_line(aes(group=rep)) +
  geom_point(aes(group=rep)) + 
  facet_wrap(~gene, nrow=1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# Compare T7Hi to deoptimized gene 10 strain
counts %>% filter(strain %in% c('T7Hi', '11-44'), gene %in% c('gene 9', 'gene 10A'), rep !='2') %>% 
  group_by(strain) %>%
  ggplot(aes(x=strain, y=tpm, fill=strain)) + 
  stat_summary(geom="bar", fun.data='mean_se') + 
  geom_line(aes(group=rep)) +
  geom_point(aes(group=rep)) + 
  facet_wrap(~gene, nrow=1) + 
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
  

# Look at relative abundance for genes 9 and 10 -- without rep2
counts %>% filter(gene %in% c('gene 9', 'gene 10A'), rep != '2') %>%
  group_by(strain) %>%
  ggplot(aes(x=background, y=tpm, fill=background)) + 
  stat_summary(geom="bar", fun.data='mean_se') + 
  geom_point(aes(group=rep), show.legend = FALSE) + 
  facet_grid(gene ~ knockout, labeller = labeller(gene=label_value, knockout=label_parsed), scales='free_x', space='free_x') + 
  scale_x_discrete(labels=back_labs) +
  panel_border() +
  labs(x='genetic background') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.text.align = 0) 






