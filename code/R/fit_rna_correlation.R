# Calculate correlation between fitness and rna abundance for each gene
library(tidyverse)
library(cowplot)
library(openxlsx)
library(broom)
library(stringr)

fit <- read.xlsx("../../data/results/Fitness.xlsx")
rna <- read.csv("../../data/results/counts_rna_abundance.csv")

rna$gene <- factor(rna$gene, levels=unique(rna$gene), ordered=TRUE)

# mean fitness for each strain 
mean_fit <- fit %>% group_by(strain_rna) %>% summarize_at(vars(Fitness), funs(mean_fit=mean)) %>% 
  dplyr::rename(strain=strain_rna, mean_fit=mean_fit)
mean_fit$strain <- factor(mean_fit$strain)

# Mean rna abundance for each gene in each strain
mean_rna <- rna %>% group_by(strain, gene) %>% summarize_at(vars(tpm), funs(mean_count=mean))
# consolidate data frames
mean_data <- left_join(mean_rna, mean_fit)

# This was used to calculate correlation using the individual RNA counts instead of the mean counts
# data <- left_join(rna, mean_fit)

res <- mean_data %>% group_by(gene) %>% 
  summarize(p.value=cor.test(mean_count, mean_fit)$p.value,
            cor=cor.test(mean_fit, mean_count)$estimate) %>% 
  mutate(R2 = round(cor*cor, 3),
         fdr = p.adjust(p.value, method = "fdr")) %>% 
  arrange(desc(R2))

res <- res %>% select(gene, R2) %>% filter(!gene=="10B")
write.csv(res, "../../data/results/fitness_rna_correlation.csv", row.names = F)

# `data` should not be here?
final <- left_join(mean_data, res) # %>% 
  # left_join(data)

# Plot RNA abundance vs. fitness and label with R-squared values
corr_plot <- final %>% filter(gene %in% c('8', '9', '10A', '11', '12')) %>%
  ggplot(aes(x=mean_fit, y=mean_count)) + 
  geom_point() + 
  geom_text(aes(x=39, y=-Inf, label=paste("R^{2} ==", R2)), vjust=-0.5, parse=T) +
  geom_smooth(method="lm") +
  scale_x_continuous(expand=c(0,0), limits=c(25.5,46.5)) +
  labs(x="fitness (doublings/hour)", y = "RNA abundance (tpm)") +
  facet_wrap(~gene, scales="free_y", nrow=1,
             labeller=labeller(gene=as_labeller(function(string, prefix="gene") paste(prefix, string))))

corr_plot

save_plot("../../figures/corr_plot.pdf", corr_plot, base_width=9, base_height=3)

