# Bootstrap 10A counts for single vs. double knockout
library(tidyverse)
library(cowplot)

rna <- read.csv("../../data/results/counts_rna_abundance.csv")
rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)


wt <- rna %>% filter(strain=="T7Hi", gene=="10A")
phi9 <- rna %>% filter(strain=='phi9v2', gene=="10A")
phi10 <- rna %>% filter(strain=='phi10', gene=="10A")
phi910 <- rna %>% filter(strain=='phi910v2', gene=="10A")

df <- data.frame(bootstrap= 1:10000)
df %>% group_by(bootstrap) %>%
  mutate(p9=(mean(sample(wt$tpm, length(wt$tpm), replace=TRUE)) - mean(sample(phi9$tpm, length(phi9$tpm), replace=TRUE))),
           p10=(mean(sample(wt$tpm, length(wt$tpm), replace=TRUE)) - mean(sample(phi10$tpm, length(phi10$tpm), replace=TRUE))), 
         double=(mean(sample(wt$tpm, length(wt$tpm), replace=TRUE)) - mean(sample(phi910$tpm, length(phi910$tpm), replace=TRUE)))) %>% 
  mutate(difference = double-(p9+p10)) %>%
  ungroup() -> bootstraps

# Calculate percentage that difference between dbl and single knockouts together is greater than 0
bootstraps %>% tally(difference>0)/length(bootstraps$difference)

grp <- data.frame(promoter=c('p9', 'p10', 'double'),
                  grp=c('single', 'single', 'double'))


bootstraps %>% select(bootstrap, p9, p10, double) %>% gather(promoter, loss, p9:double) %>% 
  group_by(promoter) %>% summarize(mean_loss=mean(loss)) %>% left_join(grp) %>%
  ggplot(aes(x=grp, y=mean_loss, fill=promoter)) + 
    geom_bar(stat='identity')

