# Differential expression analysis 
library(tidyverse)
library(cowplot)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

tpm_stats <- read.csv("../../data/results/pairwise_t_vs_wt.csv")
# deseq2_stats are not used as fare as I can tell
deseq2_stats <- read.csv("../../data/results/deseq2_results_vs_wt.csv")

rna <- read.csv("../../data/results/counts_rna_abundance.csv")

tpm_stats %>% filter(strain %in% c("8st-910evo", "910L2evo"), gene %in% c("8", "9", "10A", "11", "12"))

# take multiple samples from counts data and calculate means for expression data and convert to ratios relative to reference strain
calc_mean_ratio <- function(g, ref, treat, rna_df){
  treatment <- rna_df %>% filter(strain==treat, gene==g)
  control <- rna_df %>% filter(strain==ref, gene==g)
  df <- data.frame(bootstrap = 1:1000)
  df %>% group_by(bootstrap) %>% 
    mutate(ratio = mean(sample(treatment$tpm, length(treatment$tpm), replace=TRUE))/mean(sample(control$tpm, length(control$tpm), replace=TRUE))) %>%
    ungroup() %>%
    summarize(mean_ratio = mean(ratio), sd = sd(ratio)) %>% mutate(strain=treat, gene=g) -> bootstraps
  return(bootstraps)
}
calc_all_ratios <- function(treat.list, gene.list, rna_df, ref){
  map_df(treat.list, function(treat){
    map_df(gene.list, calc_mean_ratio, rna_df=rna_df, ref=ref, treat=treat)
  })
}

strains.list=unique(counts$strain)
# fold change in tpm compared to wt strain
ratio_wt <- calc_all_ratios(treat.list=strains.list, gene.list=c('8', '9', '10A', '11', '12'), rna_df=counts, ref='T7Hi') %>% 
  mutate(log2fc=log2(mean_ratio)) 
# fold change in tpm compared to 10deop strain
ratio_deop <- calc_all_ratios(treat.list=c('11-44', '11-44-phi9v2', '11-44-phi10'), 
                              gene.list=c('8', '9', '10A', '11', '12'), rna_df=counts, ref='11-44')


wt_ratio <- inner_join(tpm_stats %>% filter(control=='T7Hi'), ratio_wt) %>% 
  filter(!strain %in% c('11-44', '8st-910evo', '910L2evo'))
deop_ratio <- inner_join(tpm_stats %>% filter(control=='11-44'), ratio_deop)
df <- data.frame(strain='11-44-910v2',
                 background='10deop',
                 knockout='phi910',
                 gene='10A')

deop_ratio <- bind_rows(deop_ratio, df)

wt_ratio$gene <- factor(wt_ratio$gene, 
                      levels=c('8', '9', '10A', '11', '12'))
wt_ratio$background <- factor(wt_ratio$background,
                              levels=c('wt', '10deop', '8st'),
                              labels=c(
                                'wt' = 'wt',
                                '10deop' = expression(paste('10'['deop'])),
                                '8st' = expression(paste('8'[Delta]['stop']))
                              ))
wt_ratio$knockout <- factor(wt_ratio$knockout, 
                            levels=c('phi9', 'phi10', 'phi910'),
                            labels=c(
                              'phi9' = expression(paste(Delta, phi, '9')),
                              'phi10' = expression(paste(Delta, phi, '10')),
                              'phi910' = expression(paste(Delta, phi, '9/', phi, '10'))
                            ))
deop_ratio$gene <- factor(deop_ratio$gene, 
                          levels=c('8', '9', '10A', '11', '12'))
deop_ratio$background <- factor(deop_ratio$background,
                              labels=c(
                                '10deop' = expression(paste('10'['deop']))
                              ))
deop_ratio$knockout <- factor(deop_ratio$knockout,
                              levels=c('phi9', 'phi10', 'phi910'),
                              labels=c(
                                'phi9' = expression(paste(Delta, phi, '9')),
                                'phi10' = expression(paste(Delta, phi, '10')), 
                                'phi910' = expression(paste(Delta, phi, '9/', phi, '10')))
                              )

df_nas <- data.frame(x=3,
                 y=0.75,
                 label=c('','','','','','NA','','NA',''),
                 background=c('wt', 'wt', 'wt', 
                              'paste("10"["deop"])', 'paste("10"["deop"])', 'paste("10"["deop"])', 
                              'paste("8"[Delta]["stop"])', 'paste("8"[Delta]["stop"])', 'paste("8"[Delta]["stop"])'),
                 knockout=c('paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")', 
                            'paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")', 
                            'paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")'))

fc_wt <- wt_ratio %>%
  ggplot(aes(x=gene, y=mean_ratio, fill=background)) + 
  geom_hline(aes(yintercept=1), color="#E69F00") +
  # annotate("text", x=5, y=1.1, label="wt", size=3.5) + 
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymin=mean_ratio-(2*sd), ymax=mean_ratio+(2*sd)), width=0.35) +
  geom_text(data=df_nas, aes(x=x, y=y, label=label), size=10, color="grey", inherit.aes=F) + 
  facet_grid(background~knockout, labeller=label_parsed) + 
  scale_fill_manual(values=cbPalette) + 
  panel_border() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1.5)) + 
  labs(x='gene', y='relative RNA \nabundance') +
  theme(legend.position = 'none')

fc_wt

data.segm<-data.frame(x=-Inf,y=1,xend=Inf,yend=1,
                      background='paste("10"["deop"])',knockout=c('paste(Delta, phi, "9")', 'paste(Delta, phi, "10")'))
df_nas_2 <- data.frame(x=3,
                     y=0.75,
                     label=c('','','NA'),
                     background=c('paste("10"["deop"])'), 
                     knockout=c('paste(Delta, phi, "9")', 'paste(Delta, phi, "10")', 'paste(Delta, phi, "9/", phi, "10")'))
fc_deop <- deop_ratio %>%
  ggplot(aes(x=gene, y=mean_ratio)) + 
  # annotate("text", x=5, y=1.1, label="wt", size=3.5) + 
  geom_bar(stat='identity', fill="#56B4E9") + 
  geom_errorbar(aes(ymin=mean_ratio-sd, ymax=mean_ratio+sd), width=0.35) + 
  geom_segment(data=data.segm, aes(x=x, xend=xend, y=y, yend=yend), inherit.aes=F, color='#56B4E9') +
  geom_text(data=df_nas_2, aes(x=x, y=y, label=label), size=10, color="grey", inherit.aes=F) + 
  facet_grid(background~knockout, labeller=label_parsed) +
  panel_border() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1.5)) + 
  labs(x='gene', y='relative RNA \nabundance') +
  theme(legend.position = 'none')

fc_deop

plot_grid(fc_wt, fc_deop, labels="AUTO", ncol=1, align='v', axis='r', rel_heights=c(1,.53)) -> plots

plots

save_plot("../../figures/fc_plots.pdf", plots, base_width = 6, base_height = 6)
