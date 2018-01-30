# Analyzing RNA count data from Hisat2 alignments, normalizing by calculating TPM 
library(tidyverse)
library(cowplot)
library(broom)
library(colorspace)
library(ggrepel)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

rna <- read.csv("../../data/results/counts_rna_abundance.csv")
evo_pvals <- read.csv("../../data/results/pairwise_t_evo.csv")

# vectors for labels to be used for plots
knockout_labs <- c('wt' = 'wt',
                   'phi9' = expression(paste(Delta, phi, '9')),
                   'phi10' = expression(paste(Delta, phi, '10')),
                   'phi910' = expression(paste(Delta, phi, '9/', phi, '10')))
strain_labs <- c(
  'T7Hi' = 'wt',
  'phi9v2' = expression(paste(Delta, phi, '9'['wt'])),
  'phi910v2' = expression(paste(Delta, phi, '9/', phi, '10'['wt'])),
  'phi10' = expression(paste(Delta, phi, '10'['wt'])),
  '11-44' = expression(paste('10'['deop'])),
  '11-44-phi9v2' = expression(paste(Delta, phi, '9'['10'['deop']])),
  '11-44-phi10' = expression(paste(Delta, phi, '10'['10'['deop']])),
  '8st-9' = expression(paste(Delta, phi, '9'['8'[Delta]['stop']])),
  '8st-910' = expression(paste(Delta, phi, '9/', phi, '10'['8'[Delta]['stop']])),
  '910L2evo' = expression(paste(Delta, phi, '9/', phi, '10evo'['wt-L2'])),
  '8st-910evo' = expression(paste(Delta, phi, '9/', phi, '10evo'['8'[Delta]['stop']]))
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


# Set levels for knockout and background for labeling purposes. 
rna$knockout <- factor(rna$knockout, 
                       levels=c('wt', 'phi9', 'phi10', 'phi910'))
rna$knockout1 <- factor(rna$knockout, 
                        levels=c('wt', 'phi9', 'phi10', 'phi910'),
                        labels=c(
                          'wt' = 'wt',
                          'phi9' = expression(paste(Delta, phi, '9')),
                          'phi10' = expression(paste(Delta, phi, '10')),
                          'phi910' = expression(paste(Delta, phi, '9/', phi, '10'))))
rna$background <- factor(rna$background, 
                         levels=c('wt', '10deop', '8st'),
                         labels=back_labs)

#### Distribution across all genes
rna_t7 <- rna %>% filter(!strain %in% c("910L2evo", "8st-910evo")) %>%
  ggplot(aes(x=gene, y=tpm, fill=background)) +
  stat_summary(geom="bar", fun.data='mean_se') + 
  scale_fill_manual(values=cbPalette, 
                    labels=parse(text=back_labs)) + 
  scale_y_continuous(expand=c(0,0)) +
  panel_border() + 
  labs(y="RNA abundance (tpm)") + 
  geom_segment(aes(x=42.5, xend=42.5, y=0, yend=Inf), color='red') +
  geom_segment(aes(x=43.5, xend=43.5, y=0, yend=Inf), color='red') +
  facet_grid(background + knockout1~., scales='free_x', space='free_x',
             labeller=labeller(background=label_parsed, knockout1=label_parsed)) +
  #geom_rect(aes(xmin=41.5, xmax=48.5, ymin=0, ymax=Inf), fill="transparent", color='red') +
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5),
        legend.position="none",
        panel.spacing.y = unit(1, "lines"))

rna_t7
save_plot("../../figures/rna_distribution.pdf", rna_t7, base_width=10, base_height=10)


# Scatter plots of mean RNA abundance for mutants vs. wt
# calculate mean tpm for each gene grouped by strain, label gene 9 and 10A, output data frame
rna %>% filter(!gene=='10B', !grepl("evo|^11-44$", strain)) %>% select(rep, background, gene, strain, tpm) %>% group_by(strain, gene) %>% 
  summarize(mean_tpm = mean(tpm)) %>% ungroup() %>%
  spread(strain, mean_tpm) %>% mutate(label=if_else(grepl('^9$|^10A$', gene), paste(gene), '')) %>%
  gather(strain, tpm, phi9v2:`8st-910`) -> tpm_scatter

# lable information to be joined to tpm_scatter
labels <-data.frame(
  strain=c('phi9v2', 'phi910v2','phi10', '11-44', '11-44-phi9v2', '11-44-phi10', '8st-9', '8st-910', '8st-910evo', '910L2evo'),
  background=c('wt', 'wt', 'wt', '10deop', '10deop', '10deop', '8st', '8st', '8st', 'wt'),
  knockout=c('phi9', 'phi910', 'phi10', 'wt', 'phi9', 'phi10', 'phi9', 'phi910', 'phi910', 'phi910')
)

tpm_scatter <- tpm_scatter %>% left_join(labels)
tpm_scatter$background <- factor(tpm_scatter$background, 
                                 levels=c('wt', '10deop', '8st'),
                                 labels=back_labs)
tpm_scatter$knockout <- factor(tpm_scatter$knockout, 
                               levels=c('phi9', 'phi10', 'phi910'),
                               labels=c(
                                 'phi9' = expression(paste(Delta, phi, '9')),
                                 'phi10' = expression(paste(Delta, phi, '10')),
                                 'phi910' = expression(paste(Delta, phi, '9/', phi, '10'))))

# data frame used to label the 'NA' panels lacking data
df <- data.frame(x=0.07,
                 y=0.07,
                 label=c('','','','','','NA','','NA',''),
                 background=c('wt', 'wt', 'wt', 
                              'paste("10"["deop"])', 'paste("10"["deop"])', 'paste("10"["deop"])', 
                              'paste("8"[Delta]["stop"])', 'paste("8"[Delta]["stop"])', 'paste("8"[Delta]["stop"])'),
                 knockout=c('paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")', 
                            'paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")', 
                            'paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")'))

# Plot tpm_scatter by knockout and genetic background
scatter_plot <- tpm_scatter %>% filter(!grepl("evo|^11-44$", strain)) %>%
  mutate(T7Hi=T7Hi/1e6, tpm=tpm/1e6) %>%
  ggplot(aes(x=T7Hi, y=tpm)) +
  geom_segment(aes(x=-Inf, xend=Inf, y=-Inf, yend=Inf), inherit.aes=F, color="grey") + 
  geom_text_repel(aes(label = label), nudge_y = -0.005, segment.color="black") +
  geom_point(aes(color=background)) +
  geom_text(data=df, aes(x=x, y=y, label=label), size=10, color="grey", inherit.aes=F) + 
  facet_grid(background~knockout, labeller = label_parsed) +
  scale_color_manual(values=cbPalette) + 
  panel_border() + 
  labs(x="RNA abundance \n(wild-type)", y="RNA abundance \n(mutant)") + 
  xlim(0, 0.14) + 
  ylim(0, 0.14) + 
  theme(legend.position="none") 

scatter_plot

save_plot("../../figures/rna_scatter.pdf", scatter_plot, base_width = 7, base_height=7)



# Plot subset of genes (8, 9, 10A, 11, 12)
# dummy data used to allow all points to be in-frame
dummy2 <- data.frame(strain="T7Hi",
                     gene=c("8", "9", "10A", "11", "12"),
                     yval=c(0.048, 0.08, 0.165, 0.018, 0.0125))

abundance_subset_plot <- rna %>% filter(gene %in% c("8", "9", "10A", "11", "12"), !strain %in% c("910L2evo", "8st-910evo")) %>%
  ggplot(aes(x=knockout, y=tpm, fill=background)) + 
  stat_summary(geom="bar", fun.data = 'mean_se') +
  geom_blank(data=dummy2, aes(x=1, y=yval), inherit.aes=F) + 
  geom_point(show.legend = FALSE) +
  geom_line(aes(group=rep), show.legend=FALSE) +
  scale_fill_manual(values=cbPalette, 
                    labels=parse(text=back_labs)) + 
  scale_x_discrete(labels=knockout_labs) +
  scale_y_continuous(expand=c(0,0)) +
  panel_border() +
  facet_grid(gene~background, labeller = labeller(gene=as_labeller(function(string, prefix='gene') paste(prefix, string)), 
                                                  background=label_parsed), scales='free', space='free_x') + 
  labs(x = "promoter knockout", y = "RNA abundance (tpm)") + 
  theme(legend.position='none',
        legend.title = element_blank(),
        legend.key.size = unit(0.85, 'cm'),
        legend.text.align = 0,
        panel.spacing.y = unit(1, "lines"),
        axis.text.x = element_text(angle = -45, hjust=0, vjust=1))

abundance_subset_plot
save_plot("../../figures/sub_genes_plot.pdf", abundance_subset_plot, base_width=6, base_height=7)



# Initial and evolved RNA abundances for genes 8-12
# info for labelling plot
line_labs <- c(
  'T7Hi' = 'wt',
  'phi910v2' = 'initial',
  '8st-910' = 'initial',
  '910L2evo' = 'evolved',
  '8st-910evo' = 'evolved'
)
lines <- data.frame(treatment=c('wt', 'wt-910', '910L2evo', '8st-910', '8st-910evo'),
                    line=c('wt', 'wt-910', 'wt-910', '8st', '8st'))
lines$line1 <- factor(lines$line, levels=c('wt', 'wt-910', '8st'),
                     labels=c('wt', 
                              'wt-910' = expression(paste(Delta, phi, '9/', phi, '10'['wt'])), 
                              '8st' = expression(paste(Delta, phi, '9/', phi, '10'['8'[Delta]['stop']]))))

# P-values for genes 8-12, label significant values with "*"
evo <- evo_pvals %>% filter(gene %in% c("8", "9", "10A", "11", "12"))
evo$star[evo$padj > .05] <- "NS"
evo$star[evo$padj <= .05]  <- "*"
evo2 <- left_join(evo, rna) %>% group_by(gene) %>% mutate(ypos=max(tpm)) 

# Combine RNA expression data with p-values
evo_data <- left_join(rna, evo2) %>% group_by(knockout) %>% mutate(xpos=length(unique(strain))/2) %>% ungroup()
evo_data$strain <- factor(evo_data$strain, 
                        levels=c('T7Hi', 'phi910v2', '910L2evo', '8st-910', '8st-910evo'))
evo_data$gene <- factor(evo_data$gene, levels=evo_data$gene[1:60])

# dummy data for plotting all points in-frame
dummy <- data.frame(strain="T7Hi",
                    gene=c("8", "9", "10A", "11", "12"),
                    yval=c(0.04, 0.08, 0.16, 0.0175, 0.0125))

evo_data %>% filter(strain %in% c("T7Hi", "phi910v2", "910L2evo", "8st-910", "8st-910evo"), 
                    gene %in% c("8", "9", "10A", "11", "12")) %>% left_join(lines) %>%
  ggplot(aes(x=strain, y=tpm, fill=strain)) +
  stat_summary(geom="bar", fun.data='mean_se', position='dodge') +
  geom_point(position=position_jitter(w=0.1,h=0), show.legend = FALSE) +
  geom_blank(data=dummy, aes(x=1, y=yval), inherit.aes=F) + 
  scale_fill_manual(values=c('#E69F00', '#E69F00', desaturate('#E69F00', 0.7), '#009E73', desaturate('#009E73', 0.7))) +
  scale_x_discrete(labels=line_labs) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(gene~line1, scales='free', space="free_x", 
             labeller=labeller(gene=as_labeller(function(string, prefix='gene') paste(prefix, string)), 
                               line1=label_parsed)) +
  geom_text(aes(x=1.5, y=ypos*1.55, label=star), size = 6, color='black', vjust=0) + 
  geom_segment(aes(x=xpos/2, xend=xpos, y=ypos*1.4, yend=ypos*1.4), color='black') +
  labs(y='RNA abundance (tpm)') +
  panel_border() + 
  theme(legend.position="none",
        panel.spacing.y = unit(1, "lines")) -> evo_rna

evo_rna
save_plot("../../figures/evo_rna.pdf", evo_rna, base_height=7, base_width=7)
