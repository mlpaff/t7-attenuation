# Analyzing RNA count data from Hisat2 alignments, normalizing by calculating TPM 
library(tidyverse)
library(cowplot)
library(broom)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

rna <- read.csv("../../data/results/hisat2_rna_abundance.csv")

# vectors for labels to be used for plots
knockout_labs <- c('wt' = 'wt',
                   'phi9' = expression(paste(Delta, phi, '9')),
                   'phi10' = expression(paste(Delta, phi, '10')),
                   'phi910' = expression(paste(Delta, phi, '9/10')))
strain_labs <- c(
  'T7Hi' = 'wt',
  'phi9v2' = expression(paste(Delta, phi, '9'['wt'])),
  'phi910v2' = expression(paste(Delta, phi, '9/10'['wt'])),
  'phi10' = expression(paste(Delta, phi, '10'['wt'])),
  '11-44' = expression(paste('10'['deop'])),
  '11-44-phi9v2' = expression(paste(Delta, phi, '9'['10'['deop']])),
  '11-44-phi10' = expression(paste(Delta, phi, '10'['10'['deop']])),
  '8st-9' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])),
  '8st-910' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']])),
  '910L2evo' = expression(paste(Delta, phi, '9/10evo'['wt-L2'])),
  '8st-910evo' = expression(paste(Delta, phi, '9/10evo'['8'[Delta]['stop']]))
)
back_labs <- c(
  'wt' = 'wt',
  '10deop'=expression(paste('10'['deop'])),
  '8st'=expression(paste('8'[Delta]['stop']))
)
# Set levels for genes and strain
rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)
rna$strain <- factor(rna$strain, 
                     levels=c('T7Hi', 'phi9v2', 'phi10', 'phi910v2', '910L2evo', 
                              '11-44', '11-44-phi9v2', '11-44-phi10', 
                              '8st-9', '8st-910', '8st-910evo'))

# Normalize the raw RNA counts, first to the sum of raw counts for all genes up to and including gene 7.7, then by TPM
rna %>% group_by(rep, strain) %>% 
  mutate(nfactor=sum(counts[gene<='7.7'])) %>%
  mutate(normcounts=counts/nfactor) %>%
  mutate(rpk=normcounts/((stop-start)/1000)) %>%
  mutate(rpm=sum(rpk)) %>%
  mutate(tpm=rpk/rpm) -> counts

counts$gene <- factor(counts$gene, levels=rna$gene[1:60])


# Set levels for knockout and background for labeling purposes. 
counts$knockout <- factor(counts$knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'))
counts$knockout1 <- factor(counts$knockout, 
                        levels=c('wt', 'phi9', 'phi10', 'phi910'),
                        labels=c(
                          'wt' = 'wt',
                          'phi9' = expression(paste(Delta, phi, '9')),
                          'phi10' = expression(paste(Delta, phi, '10')),
                          'phi910' = expression(paste(Delta, phi, '9/10'))))
counts$background <- factor(counts$background, 
                         levels=c('wt', '10deop', '8st'),
                         labels=back_labs)

# Take a look at the distributions for rep 1 and rep 2   
counts %>% filter(rep %in% c('2', '3')) %>%
  ggplot(aes(x=gene, y=tpm, fill=background)) +
  stat_summary(geom="bar", fun.data='mean_se') +
  geom_point(aes(group=rep, color=rep)) + 
  scale_fill_manual(values=cbPalette, 
                    labels=parse(text=back_labs)) + 
  scale_y_continuous(expand=c(0,0)) +
  panel_border() + 
  labs(y="relative mRNA abundance") + 
  facet_grid(strain~., scales='free_x', space='free_x') +
  #geom_rect(aes(xmin=41.5, xmax=48.5, ymin=0, ymax=Inf), fill="transparent", color='red') +
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5),
        legend.position="none",
        panel.spacing.y = unit(1, "lines"))

#### Distribution across all genes
rna_t7 <- counts %>% filter(!strain %in% c("910L2evo", "8st-910evo")) %>%
  ggplot(aes(x=gene, y=tpm, fill=background)) +
  stat_summary(geom="bar", fun.data='mean_se') + 
  scale_fill_manual(values=cbPalette, 
                    labels=parse(text=back_labs)) + 
  scale_y_continuous(expand=c(0,0)) +
  panel_border() + 
  labs(y="relative mRNA abundance") + 
  facet_grid(background + knockout1~., scales='free_x', space='free_x',
             labeller=labeller(background=label_parsed, knockout1=label_parsed)) +
  #geom_rect(aes(xmin=41.5, xmax=48.5, ymin=0, ymax=Inf), fill="transparent", color='red') +
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5),
        legend.position="none",
        panel.spacing.y = unit(1, "lines"))

rna_t7
save_plot("../../figures/rna_distribution.pdf", rna_t7, base_width=10, base_height=10)


# Distribution for all genes for initial and evolved lines
line_labs <- c(
  'T7Hi' = 'wt',
  'phi910v2' = 'initial',
  '8st-910' = 'initial',
  '910L2evo' = 'evolved',
  '8st-910evo' = 'evolved'
)
lines <- data.frame(treatment=c('wt', 'wt-910', '910L2evo', '8st-910', '8st-910evo'),
                    line=c('wt', 'wt-910', 'wt-910', '8st', '8st'))
lines$line <- factor(lines$line, levels=c('wt', 'wt-910', '8st'), 
                     labels=c('wt', 
                              'wt-910' = expression(paste(Delta, phi, '9/10'['wt'])), 
                              '8st' = expression(paste(Delta, phi, '9/10'['8'[Delta]['stop']]))))

counts %>% filter(strain %in% c("T7Hi", "phi910v2", "910L2evo", "8st-910", "8st-910evo"), 
                  gene %in% c("8", "9", "10A", "11", "12", "13")) %>% left_join(lines) %>%
  ggplot(aes(x=strain, y=tpm, fill=line)) +
  stat_summary(geom="bar", fun.data='mean_se', position='dodge') +
  #geom_jitter(width = 0.1, height = 0) + 
  geom_point(position=position_jitter(w=0.1,h=0), show.legend = FALSE) +
  scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=line_labs) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(gene~line, scales='free', space="free_x", 
             labeller=labeller(gene=as_labeller(function(string, prefix='gene') paste(prefix, string)), 
                               line=label_parsed)) +
  labs(y='relative mRNA abundance') +
  panel_border() + 
  theme(legend.position="none",
        panel.spacing.y = unit(1, "lines")) -> evo_rna

evo_rna
save_plot("../../figures/evo_rna.pdf", evo_rna, base_height=7, base_width=7)


counts %>% filter(strain %in% c("T7Hi", "phi910v2", "910L2evo", "8st-910", "8st-910evo")) %>% left_join(lines) %>%
  ggplot(aes(x=gene, y=tpm, fill=line)) +
  stat_summary(geom="bar", fun.data='mean_se', position='dodge') +
  #geom_jitter(width = 0.1, height = 0) + 
  #geom_point(position=position_jitter(w=0.1,h=0), show.legend = FALSE) +
  scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=line_labs) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(strain~., scales='free_x', space="free") +
  labs(y='relative mRNA abundance') +
  panel_border() + 
  theme(legend.position="none",
        panel.spacing.y = unit(1, "lines"),
        axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5))


# Look at a subset of genes (8, 9, 10A, 11, 12, 13)
abundance_subset_plot <- counts %>% filter(gene %in% c("8", "9", "10A", "11", "12", "13"), !strain %in% c("910L2evo", "8st-910evo")) %>%
  ggplot(aes(x=knockout, y=tpm, fill=background)) + 
  stat_summary(geom="bar", fun.data = 'mean_se') +
  geom_point(show.legend = FALSE) +
  geom_line(aes(group=rep), show.legend=FALSE) +
  scale_fill_manual(values=cbPalette, 
                    labels=parse(text=back_labs)) + 
  scale_x_discrete(labels=knockout_labs) +
  scale_y_continuous(expand=c(0,0)) +
  panel_border() +
  facet_grid(gene~background, labeller = labeller(gene=as_labeller(function(string, prefix='gene') paste(prefix, string)), 
                                                  background=label_parsed), scales='free', space='free_x') + 
  labs(x = "promoter knockout", y = "relative mRNA abundance") + 
  theme(legend.position='none',
    legend.title = element_blank(),
    legend.key.size = unit(0.85, 'cm'),
    legend.text.align = 0,
    panel.spacing.y = unit(1, "lines"),
    axis.text.x = element_text(angle = -45, hjust=0, vjust=1))

abundance_subset_plot
save_plot("../../figures/sub_genes_plot.pdf", abundance_subset_plot, base_width=6, base_height=7)
